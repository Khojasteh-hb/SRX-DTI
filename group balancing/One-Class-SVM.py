# Data balancing with One-Class SVM

import numpy as np
import pandas as pd
from sklearn.svm import OneClassSVM

import warnings
warnings.filterwarnings('ignore')

# *** groups ***
# Read drug-target interactions
df_inter = pd.read_csv('df_NR_AB.csv')


df_inter.drop(df_inter.columns[df_inter.columns.str.contains('unnamed', case=False)], axis=1, inplace=True)
#print(df_inter.head())

count_label_0 = sum(df_inter['label'] == 0)
count_label_1 = sum(df_inter['label'] == 1)


# **Splitting samples of label(0) **
df_label_0 = df_inter.loc[df_inter['label'] == 0]
df_label_1 = df_inter.loc[df_inter['label'] == 1]

df_X = df_label_0.copy()
df_X.drop(['drug_no', 'protein_no', 'label'], axis=1, inplace=True)

print(df_X.shape[0])

svm = OneClassSVM(kernel='rbf', gamma=1.0/df_X.shape[0], tol=0.001, nu=0.5, shrinking=True, cache_size=80)
svm = svm.fit(df_X.values)
print(svm)

# decision_function it says that its the distance between the hyperplane and the test instance
# flatten() function we can flatten a matrix to one dimension in python.
scores = svm.decision_function(df_X.values).flatten()
maxvalue = np.max(scores)
print(maxvalue)

scores = maxvalue - scores
print(scores)

scores = pd.DataFrame(scores)
scores.rename(columns={0: 'score'}, inplace=True)

samples = [i for i in range(count_label_0)]
samples = pd.DataFrame(samples)
samples.rename(columns={0: 'sample'}, inplace=True)

df_scores = pd.concat([samples, scores], axis=1, sort=False)
df_scores_reordered = df_scores.sort_values('score', ascending=True)

df_low_score = df_scores_reordered.head(count_label_1)
print(df_low_score.head())
print(len(df_low_score))

index_list = df_low_score['sample'].tolist()
print(index_list)

df_final_label_0 = df_label_0.iloc[index_list, :]
print(df_final_label_0.head())

# append method
df_final = df_final_label_0.append(df_label_1)
df_final.reset_index(inplace=True)
df_final.drop(['index'], axis=1, inplace=True)
print(df_final.head())
print(len(df_final))


# *** groups ***
df_final.to_csv('df_final_NR_AB.csv')



