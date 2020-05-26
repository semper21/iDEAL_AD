'''
Created on May 11, 2020

@author: ywkim
'''
import numpy as np
import pandas as pd
from pathlib import Path

from germline_analyses import get_matrix_subset, get_list_from_csv

from sklearn import metrics
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import AdaBoostClassifier

from sklearn.svm import SVC

if __name__ == '__main__':
    output_folder =  str(Path().absolute()) + '/output/'
    validation_folder = str(Path().absolute()) + '/output_ADSP_extension/'

    gene_file = output_folder + 'iDEAL_genelist.txt'
    gene_list = get_list_from_csv(gene_file, 'Gene', sep='\t')

    # these are gene x patient matrices
    train_ADe2 = get_matrix_subset(output_folder, 'sum_ea_matrix_ADe2.tsv', gene_list, sep='\t', index_col=0)
    train_HCe4 = get_matrix_subset(output_folder, 'sum_ea_matrix_HCe4.tsv', gene_list, sep='\t', index_col=0)
    val_ADe2 = get_matrix_subset(validation_folder, 'sum_ea_matrix_ADe2.tsv', gene_list, sep='\t', index_col=0)
    val_HCe4 = get_matrix_subset(validation_folder, 'sum_ea_matrix_HCe4.tsv', gene_list, sep='\t', index_col=0)

    x = np.vstack((train_ADe2.T, train_HCe4.T))
    y = [1] * len(train_ADe2.T) + [0] * len(train_HCe4.T)

    x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.25, random_state=0)

    """DecisionTree as base_estimator"""
    """
    max_depths = [2, 5, 10]
    learning_rates = [0.01, 0.1, 0.5]
    n_estimators = [50, 100, 200, 500]
    acc_table = pd.DataFrame(columns=['max_depth', 'learning rate', 'n_estimators', 'Train_acc', 'Train_prec',
                                      'Train_recall', 'f1', 'Train_conf_mat'], index=range(36))
    j=0
    for d in max_depths:
        for lr in learning_rates:
            for n in n_estimators:
                classifier = AdaBoostClassifier(DecisionTreeClassifier(max_depth=d), n_estimators=n,
                                                learning_rate=lr, random_state=0)
                classifier.fit(x_train, y_train)
                y_pred = classifier.predict(x_test)
                cfm = confusion_matrix(y_test, y_pred, labels=[1, 0])

                p = metrics.precision_score(y_test, y_pred)
                r = metrics.recall_score(y_test, y_pred)
                f = 2 * (p * r)/(p + r)

                acc_table.iloc[j, 0] = d
                acc_table.iloc[j, 1] = lr
                acc_table.iloc[j, 2] = n
                acc_table.iloc[j, 3] = metrics.accuracy_score(y_test, y_pred)
                acc_table.iloc[j, 4] = p
                acc_table.iloc[j, 5] = r
                acc_table.iloc[j, 6] = f
                acc_table.iloc[j, 7] = str(cfm)
                j+=1

                # print(d, lr, n, metrics.accuracy_score(y_val, y_pred), p, r, f)

    acc_table.to_csv(output_folder + 'adaboost_dtree_test_accuracy.tsv', sep='\t', index=False)
    """

    """SVC as base_estimator"""

    """
    learning_rates = [0.01, 0.1, 0.5]
    n_estimators = [100, 200, 300]
    c_range = [0.1, 1, 10, 100]
    acc_table = pd.DataFrame(columns=['learning rate', 'n_estimators', 'Train_acc', 'Train_prec',
                                      'Train_recall', 'f1', 'Train_conf_mat'], index=range(36))
    j=0
    for lr in learning_rates:
        for n in n_estimators:
            for c in c_range:
                svc = SVC(probability=True, gamma='auto', C=c, kernel='rbf')
                classifier = AdaBoostClassifier(base_estimator=svc, n_estimators=n,
                                                learning_rate=lr, random_state=0)
                classifier.fit(x_train, y_train)
                y_pred = classifier.predict(x_test)
                cfm = confusion_matrix(y_test, y_pred, labels=[1, 0])

                p = metrics.precision_score(y_test, y_pred)
                r = metrics.recall_score(y_test, y_pred)
                f = 2 * (p * r)/(p + r)

                acc_table.iloc[j, 0] = lr
                acc_table.iloc[j, 1] = n
                acc_table.iloc[j, 2] = metrics.accuracy_score(y_test, y_pred)
                acc_table.iloc[j, 3] = p
                acc_table.iloc[j, 4] = r
                acc_table.iloc[j, 5] = f
                acc_table.iloc[j, 6] = str(cfm)
                j+=1

                print(c, lr, n, metrics.accuracy_score(y_test, y_pred), p, r, f)

    acc_table.to_csv(output_folder + 'adaboost_rbf_svc_test_accuracy.tsv', sep='\t', index=False)

    """
    x_val = np.vstack((val_ADe2.T, val_HCe4.T))
    y_val = [1] * len(val_ADe2.T) + [0] * len(val_HCe4.T)

    svc = SVC(probability=True, kernel='linear')
    classifier = AdaBoostClassifier(base_estimator=svc, n_estimators=200,
                                    learning_rate=0.01, random_state=0)
    classifier.fit(x, y)
    y_pred_val = classifier.predict(x_val)
    cfm2 = confusion_matrix(y_val, y_pred_val, labels=[1, 0])

    print('validation confusion matrix = ', cfm2)
    print("Accuracy:", metrics.accuracy_score(y_val, y_pred_val))
    print("Precision:", metrics.precision_score(y_val, y_pred_val))
    print("Recall:", metrics.recall_score(y_val, y_pred_val))
