
import numpy as np
import pandas as pd
from sklearn.svm import OneClassSVM

from sklearn.model_selection import train_test_split

from sklearn import metrics

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import TimeSeriesSplit, GridSearchCV, RandomizedSearchCV

import warnings
warnings.filterwarnings('ignore')

def One_SVM_US(df):

    count_label_0 = sum(df['label'] == 0)
    count_label_1 = sum(df['label'] == 1)

    # **Splitting samples of label(0) **
    df_label_0 = df.loc[df['label'] == 0]
    df_label_1 = df.loc[df['label'] == 1]

    df_X = df_label_0.copy()
    df_X.drop(['drug_no', 'protein_no', 'label'], axis=1, inplace=True)

    print(df_X.shape[0])

    svm = OneClassSVM(kernel='rbf', gamma=1.0 / df_X.shape[0], tol=0.001, nu=0.5, shrinking=True, cache_size=80)
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

    return df_final

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

def FFS_RF(df_group):
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


def main():
    # *** groups ***
    # Read drug-target interactions
    df_inter = pd.read_csv('df_NR_AB.csv')

    df_inter.drop(df_inter.columns[df_inter.columns.str.contains('unnamed', case=False)], axis=1, inplace=True)
    # print(df_inter.head())
    df_final = One_SVM_US(df_inter)

    FFS_RF(df_final)


if __name__ == "__main__":
    main()