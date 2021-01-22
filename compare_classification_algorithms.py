'''
Created on May 8, 2020

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
from sklearn.model_selection import train_test_split

from sklearn.neural_network import MLPClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis



if __name__ == '__main__':
    output_folder =  str(Path().absolute()) + '/output_ADSP_discovery/'

    gene_file = output_folder + 'iDEAL_genelist.txt'
    gene_list = get_list_from_csv(gene_file, 'Gene', sep='\t')

    # these are gene x patient matrices
    train_ADe2 = get_matrix_subset(output_folder, 'sum_ea_matrix_ADe2.tsv', gene_list, sep='\t', index_col=0)
    train_HCe4 = get_matrix_subset(output_folder, 'sum_ea_matrix_HCe4.tsv', gene_list, sep='\t', index_col=0)
    val_ADe3 = get_matrix_subset(output_folder, 'sum_ea_matrix_ADe3.tsv', gene_list, sep='\t', index_col=0)
    val_HCe3 = get_matrix_subset(output_folder, 'sum_ea_matrix_HCe3.tsv', gene_list, sep='\t', index_col=0)

    x = np.vstack((train_ADe2.T, train_HCe4.T))
    y = [1] * len(train_ADe2.T) + [0] * len(train_HCe4.T)
    print(x.shape)

    x_val = np.vstack((val_ADe3.T, val_HCe3.T))
    y_val = [1] * len(val_ADe3.T) + [0] * len(val_HCe3.T)

    """L1-based feature selection"""
    from sklearn.feature_selection import SelectFromModel
    from sklearn.linear_model import LogisticRegression
    lr = LogisticRegression(penalty='l1', C=100, solver='liblinear', random_state=0).fit(x,y)
    model = SelectFromModel(lr, prefit=True)
    x_new = model.transform(x)
    print(x_new.shape)

    feature_idx = model.get_support(indices=True)
    x_val_new = x_val[:,feature_idx]

    names = ["Nearest Neighbors", "Linear SVM", "RBF SVM", "Gaussian Process",
             "Decision Tree", "Random Forest", "Neural Net", "AdaBoost",
             "Naive Bayes", "QDA"]

    classifiers = [
        KNeighborsClassifier(3),
        SVC(kernel="linear", C=0.025),
        SVC(gamma=2, C=1),
        GaussianProcessClassifier(1.0 * RBF(1.0)),
        DecisionTreeClassifier(),
        RandomForestClassifier(),
        MLPClassifier(alpha=1, max_iter=1000),
        AdaBoostClassifier(),
        GaussianNB(),
        QuadraticDiscriminantAnalysis()]

    acc_table = pd.DataFrame(columns=['Algorithm', 'Accuracy', 'Precision', 'Recall', 'f1', 'Confusion_matrix'])
    acc_table['Algorithm'] = names

    j = 0
    for name, classifier in zip(names, classifiers):
        classifier.fit(x_new, y)
        # y_pred_val = classifier.predict(x_val)
        y_pred_val = classifier.predict(x_val_new) # only using features selected prior
        cfm = confusion_matrix(y_val, y_pred_val, labels=[1,0])

        p = metrics.precision_score(y_val, y_pred_val)
        r = metrics.recall_score(y_val, y_pred_val)
        acc_table.iloc[j, 1] = metrics.accuracy_score(y_val, y_pred_val)
        acc_table.iloc[j, 2] = p
        acc_table.iloc[j, 3] = r
        acc_table.iloc[j, 4] = 2 * (p * r)/(p + r)
        acc_table.iloc[j, 5] = str(cfm)
        j+=1

        # plt.clf()
        # sns.heatmap(cfm, annot=True, fmt="g", cmap='Spectral_r')
        # plt.savefig(output_folder + name + '_val_apoe3_confusion_matrix.png', dpi=200)
        print(metrics.accuracy_score(y_val, y_pred_val))
    acc_table.to_csv(output_folder + 'l1fs_classification_accuracy.tsv', sep='\t', index=False)

