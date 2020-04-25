"""
Created on Apr 25, 2020

@author: ywkim


"""


import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt

from germline_analyses import get_matrix_subset

from sklearn import metrics
from sklearn.metrics import confusion_matrix
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split


if __name__ == '__main__':
    discovery_folder = str(Path().absolute()) + '/output/'
    validation_folder =  str(Path().absolute()) + '/output_ADSP_extension/'

    gene_file = discovery_folder + 'iDEAL_genelist.txt'
    df_gene = pd.read_csv(gene_file, index_col=None)
    gene_list = df_gene['Gene'].values.tolist()

    train_ADe2 = get_matrix_subset(discovery_folder, 'sum_ea_matrix_ADe2.tsv', gene_list, sep='\t', index_col=0)
    train_HCe4 = get_matrix_subset(discovery_folder, 'sum_ea_matrix_HCe4.tsv', gene_list, sep='\t', index_col=0)
    test_ADe2 = get_matrix_subset(validation_folder, 'sum_ea_matrix_ADe2.tsv', gene_list, sep='\t', index_col=0)
    test_HCe4 = get_matrix_subset(validation_folder, 'sum_ea_matrix_HCe4.tsv', gene_list, sep='\t', index_col=0)


    """Logistic Regression"""
    x = np.vstack((train_ADe2.T, train_HCe4.T))
    y = [1] * len(train_ADe2.T) + [0] * len(train_HCe4.T)

    x_train, x_test, y_train, y_test = train_test_split(x, y, test_size = 0.2)

    lr = LogisticRegression(random_state=1)
    lr.fit(x_train, y_train)

    y_pred = lr.predict(x_test)
    cfm = confusion_matrix(y_test, y_pred)
    print('test confusion matrix = ', cfm)

    print("Accuracy:", metrics.accuracy_score(y_test, y_pred))
    print("Precision:", metrics.precision_score(y_test, y_pred))
    print("Recall:", metrics.recall_score(y_test, y_pred))

    x_val = np.vstack((test_ADe2.T, test_HCe4.T))
    y_val = [1] * len(test_ADe2.T) + [0] * len(test_HCe4.T)

    y_pred_val = lr.predict(x_val)
    cfm2 = confusion_matrix(y_val, y_pred_val)
    print('validation confusion matrix = ', cfm2)

    plt.clf()
    sns.heatmap(cfm2, annot=True, cmap='Spectral_r')

    print("Accuracy:", metrics.accuracy_score(y_val, y_pred_val))
    print("Precision:", metrics.precision_score(y_val, y_pred_val))
    print("Recall:", metrics.recall_score(y_val, y_pred_val))

    plt.show()