# Classifiers (SVM, RF, MLP, XGBoost) for drug-target prediction (5-fold)

import pandas as pd
import numpy as np

from xgboost import XGBClassifier
from sklearn.model_selection import KFold
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
from itertools import cycle


# Balanced with One-SVM-US ############################################################################
df_final = pd.read_csv('df_final_NR_AB.csv')

#######################################################################################################

# Remove unnamed columns
df_final.drop(df_final.columns[df_final.columns.str.contains('unnamed', case=False)], axis=1, inplace=True)

y_final = df_final['label'].values
df_final.drop(['label', 'protein_no', 'drug_no'], axis=1, inplace=True)
#print(df_final.head())

#######################################################################################################
# Read selected feature list
df_features_final = pd.read_csv('subset_accuracy_NR_AB.csv')

#######################################################################################################

features_final = df_features_final.loc[0, 'Feature']
features_final = features_final.replace(" ", "")
features_list1 = list(features_final.split(','))
print(features_list1)
print(len(features_list1))

X_final = df_final[[c for c in df_final.columns if c in features_list1]].values
#######################################################################################################

# Balanced with random undersampling ##################################################################
df_rand = pd.read_csv('df_NR_AB_rand.csv')

#######################################################################################################

# Remove unnamed columns
df_rand.drop(df_rand.columns[df_rand.columns.str.contains('unnamed', case=False)], axis=1, inplace=True)

y_rand = df_rand['label'].values
df_rand.drop(['label', 'protein_no', 'drug_no'], axis=1, inplace=True)
#print(df_rand.head())

#######################################################################################################
# Read selected feature list
df_features_rand = pd.read_csv('subset_accuracy_NR_AB_rand.csv')

#######################################################################################################

features_rand = df_features_rand.loc[0, 'Feature']
features_rand = features_rand.replace(" ", "")
features_list2 = list(features_rand.split(','))
print(features_list2)
print(len(features_list2))

X_rand = df_rand[[c for c in df_rand.columns if c in features_list2]].values
#######################################################################################################

# prepare the cross-validation procedure
cv = KFold(n_splits=5, random_state=1, shuffle=True)


# Set up plotting area
plt.figure()
lw = 2

colors = cycle(['cyan', 'indigo', 'seagreen', 'yellow', 'blue', 'darkorange'])

## Create a XGBoost Classifier##########################################################################
classifier_xgb = XGBClassifier(eval_metric='logloss', use_label_encoder=False)

## Training XGBoost Classifier with 1SVM balanced ######################################################

mean_tpr = 0.0
mean_fpr = np.linspace(0, 1, 100)
i = 0
for (train, test), color in zip(cv.split(X_final, y_final), colors):
    probas_ = classifier_xgb.fit(X_final[train], y_final[train]).predict_proba(X_final[test])
    # Compute ROC curve and area the curve
    fpr, tpr, thresholds = roc_curve(y_final[test], probas_[:, 1])
    mean_tpr += np.interp(mean_fpr, fpr, tpr)
    mean_tpr[0] = 0.0
    roc_auc = auc(fpr, tpr)

    i += 1

mean_tpr /= cv.get_n_splits(X_final, y_final)
mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)
mean_auc = format(round(mean_auc, 2), ".2f")
plt.plot(mean_fpr, mean_tpr, label="With One-SVM-US, AUROC="+str(mean_auc))

## Training XGBoost Classifier with Random balanced #####################################################

mean_tpr = 0.0
mean_fpr = np.linspace(0, 1, 100)
i = 0
for (train, test), color in zip(cv.split(X_rand, y_rand), colors):
    probas_ = classifier_xgb.fit(X_rand[train], y_rand[train]).predict_proba(X_rand[test])
    # Compute ROC curve and area the curve
    fpr, tpr, thresholds = roc_curve(y_rand[test], probas_[:, 1])
    mean_tpr += np.interp(mean_fpr, fpr, tpr)
    mean_tpr[0] = 0.0
    roc_auc = auc(fpr, tpr)

    i += 1

mean_tpr /= cv.get_n_splits(X_rand, y_rand)
mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)
mean_auc = format(round(mean_auc, 2), ".2f")
plt.plot(mean_fpr, mean_tpr, label="With Random sampling, AUROC="+str(mean_auc))


#######################################################################################################
plt.plot([0, 1], [0, 1], color="black", lw=lw, linestyle="--")
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")

#plt.title("Group AB - EN")
#plt.title("Group AB - GPCR")
#plt.title("Group AB - IC")
plt.title("Group AB - NR")

plt.legend(loc="lower right")
plt.show()