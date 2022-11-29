# FFS-RF algorithm
# Feature selection by random forest

import numpy as np
import pandas as pd

from sklearn.model_selection import train_test_split

from sklearn import metrics

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import TimeSeriesSplit, GridSearchCV, RandomizedSearchCV

import warnings
warnings.filterwarnings('ignore')

# **** Function for removing features with correlation higher than 0.8 ****
def remove_correlated_features(df):
    correlated_features = set()
    correlation_matrix = df.corr()
    for i in range(len(correlation_matrix.columns)):
        for j in range(i):
            if abs(correlation_matrix.iloc[i, j]) > 0.8:
                colname = correlation_matrix.columns[i]
                correlated_features.add(colname)
    return correlated_features

# Run one of datasets
df_group = pd.read_csv('df_final_NR_AB.csv')

# Remove unnamed columns
df_group.rename({"Unnamed: 0": "redundant"}, axis="columns", inplace=True)
df_group.drop(["redundant"], axis=1, inplace=True)
#print(df_group.head())

y = df_group['label']
df_group.drop(['label', 'protein_no', 'drug_no'], axis=1, inplace=True)
print(df_group.head())

correlated_features = remove_correlated_features(df_group)
df_group.drop(labels=correlated_features, axis=1, inplace=True)
print(df_group.head())
X = df_group

X_train, X_test, y_train, y_test = train_test_split(df_group, y, test_size=0.2, random_state=42)

n_samples, n_features = X_train.shape

n_estimators = [5, 10, 50, 100, 150, 200, 250, 300]
max_depth = [5, 10, 25, 50, 75, 100]
min_samples_leaf = [1, 2, 4, 8, 10]
min_samples_split = [2, 4, 6, 8, 10]
max_features = ["sqrt", "log2", None]

hyperparameter = {'n_estimators': n_estimators,
                  'max_depth': max_depth,
                  'min_samples_leaf': min_samples_leaf,
                  'min_samples_split': min_samples_split,
                  'max_features': max_features,
                  }

base_model_rf = RandomForestClassifier(criterion="gini", random_state=42)
n_iter_search = 30
scoring = "accuracy"
n_selected_features = 10

# selected feature set, initialized to be empty
F = []
count = 0
ddict = {}
all_F = []
all_c = []
all_acc = []
all_model = []

while count < n_selected_features:
    max_acc = 0
    for i in X_train.columns:
        if i not in F:
            F.append(i)
            X_train_tmp = X_train[F]
            acc = 0
            rsearch_cv = RandomizedSearchCV(estimator=base_model_rf,
                                            random_state=42,
                                            param_distributions=hyperparameter,
                                            n_iter=n_iter_search,
                                            cv=5,
                                            scoring=scoring,
                                            n_jobs=-1)
            rsearch_cv.fit(X_train_tmp, y_train)
            best_estimator = rsearch_cv.best_estimator_
            y_pred = best_estimator.predict(X_test[F])
            acc = metrics.accuracy_score(y_test, y_pred)
            F.pop()
            if acc > max_acc:
                max_acc = acc
                idx = i
                best_model = best_estimator

    F.append(idx)
    count += 1

    print("The current number of features: {} - Accuracy: {}%".format(count, round(max_acc * 100, 2)))

    all_F.append(np.array(F))
    all_c.append(count)
    all_acc.append(max_acc)
    all_model.append(best_model)


c = pd.DataFrame(all_c)
a = pd.DataFrame(all_acc)
f = pd.DataFrame(all_F)
f["All"] = f[f.columns[0:]].apply(
    lambda x: ', '.join(x.dropna().astype(str)), axis=1)

all_info = pd.concat([c, a, f["All"]], axis=1)
all_info.columns = ['Num_feature', 'Accuracy', 'Feature']
all_info = all_info.sort_values(by='Accuracy', ascending=False).reset_index(drop=True)


all_info.to_csv("subset_accuracy_NR_AB.csv", index=False)

f.to_csv("feature_subset_NR_AB.csv")

# Fetching the best hyperparameters
print(rsearch_cv.best_params_)

with open('best_params_NR_AB.txt', 'w+') as file:
    file.write(str(rsearch_cv.best_params_))
