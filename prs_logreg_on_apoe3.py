'''
Created on May 7, 2020

@author: ywkim
'''

import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt

from germline_analyses import get_matrix_subset, get_list_from_csv

from sklearn import metrics
from sklearn.metrics import confusion_matrix
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split


if __name__ == '__main__':
    output_folder =  str(Path().absolute()) + '/output/'

    gene_file = output_folder + 'iDEAL_genelist.txt'
    gene_list = get_list_from_csv(gene_file, 'Gene', sep='\t')

    train_ADe2 = get_matrix_subset(output_folder, 'sum_ea_matrix_ADe2.tsv', gene_list, sep='\t', index_col=0)
    train_HCe4 = get_matrix_subset(output_folder, 'sum_ea_matrix_HCe4.tsv', gene_list, sep='\t', index_col=0)
    val_ADe3 = get_matrix_subset(output_folder, 'sum_ea_matrix_ADe3.tsv', gene_list, sep='\t', index_col=0)
    val_HCe3 = get_matrix_subset(output_folder, 'sum_ea_matrix_HCe3.tsv', gene_list, sep='\t', index_col=0)

    """Logistic Regression"""
    x = np.vstack((train_ADe2.T, train_HCe4.T))
    y = [1] * len(train_ADe2.T) + [0] * len(train_HCe4.T)

    x_train, x_test, y_train, y_test = train_test_split(x, y, test_size = 0.25, random_state=0)

    acc_table = pd.DataFrame(columns=['C_parameter', 'Train_acc', 'Train_prec', 'Train_recall', 'Train_conf_mat'])

    c_param_range = [0.001, 0.01, 0.1, 1, 10, 100, 200]
    acc_table['C_parameter'] = c_param_range
    j = 0
    for i in c_param_range:
        lr = LogisticRegression(penalty='l1', C=i, solver='liblinear', random_state=0)
        lr.fit(x_train, y_train)

        y_pred = lr.predict(x_test)
        cfm = confusion_matrix(y_test, y_pred, labels=[1,0])

        acc_table.iloc[j, 1] = metrics.accuracy_score(y_test, y_pred)
        acc_table.iloc[j, 2] = metrics.precision_score(y_test, y_pred)
        acc_table.iloc[j, 3] = metrics.recall_score(y_test, y_pred)
        acc_table.iloc[j, 4] = str(cfm)
        j+=1

        print(i, metrics.accuracy_score(y_test, y_pred))
    acc_table.to_csv(output_folder + 'sum_ea_test_accuracy.tsv', sep='\t', index=False)

    x_val = np.vstack((val_ADe3.T, val_HCe3.T))
    y_val = [1] * len(val_ADe3.T) + [0] * len(val_HCe3.T)

    lr = LogisticRegression(penalty='l1', C=100, solver='liblinear', random_state=0)
    lr.fit(x, y)

    y_pred_val = lr.predict(x_val)
    cfm2 = confusion_matrix(y_val, y_pred_val, labels=[1,0])
    from IPython import embed
    # embed()
    print('validation confusion matrix = ', cfm2)

    plt.clf()
    sns.heatmap(cfm2, annot=True, fmt="g", cmap='Spectral_r')

    print("Accuracy:", metrics.accuracy_score(y_val, y_pred_val))
    print("Precision:", metrics.precision_score(y_val, y_pred_val))
    print("Recall:", metrics.recall_score(y_val, y_pred_val))

    plt.savefig(output_folder + 'val_apoe3_confusion_matrix.png', dpi=200)