'''
Created on May 11, 2020

@author: ywkim
'''
import numpy as np
import pandas as pd
from pathlib import Path

from src.germline_analyses import get_matrix_subset, get_list_from_csv

from sklearn.ensemble import AdaBoostClassifier
from sklearn.svm import SVC

from sklearn.metrics import auc
from sklearn.metrics import plot_roc_curve
from sklearn.model_selection import StratifiedKFold
from sklearn.inspection import permutation_importance

import matplotlib.pyplot as plt

if __name__ == '__main__':
    output_folder =  str(Path().absolute()) + '/output_ADSP_discovery/'

    gene_file = output_folder + 'iDEAL_genelist.txt'
    comparison = 'APOE2_AD_v_APOE4_HC'

    # gene_file = output_folder + 'pathogenic.txt'
    # comparison = 'APOE2_ADvHC_pathogenic'

    # gene_file = output_folder + 'protective.txt'
    # comparison = 'APOE4_HCvAD_protective'

    gene_list = get_list_from_csv(gene_file, 'Gene', sep='\t')

    # these are gene x patient matrices
    ADe2 = get_matrix_subset(output_folder, 'sum_ea_matrix_ADe2.tsv', gene_list, sep='\t', index_col=0)
    ADe4 = get_matrix_subset(output_folder, 'sum_ea_matrix_ADe4.tsv', gene_list, sep='\t', index_col=0)
    HCe4 = get_matrix_subset(output_folder, 'sum_ea_matrix_HCe4.tsv', gene_list, sep='\t', index_col=0)
    HCe2 = get_matrix_subset(output_folder, 'sum_ea_matrix_HCe2.tsv', gene_list, sep='\t', index_col=0)

    x = np.vstack((ADe2.T, HCe4.T))
    y = np.asarray([1] * len(ADe2.T) + [0] * len(HCe4.T))

    n_samples, n_features = x.shape

    cv = StratifiedKFold(n_splits=5)

    svc = SVC(probability=True, kernel='linear')

    # hyperparameter search
    learning_rates = [0.01, 0.1]
    n_estimators = [100, 200]

    for lr in learning_rates:
        for n in n_estimators:
            classifier = AdaBoostClassifier(base_estimator=svc, n_estimators=n,
                                            learning_rate=lr, random_state=0)
            tprs = []
            aucs = []
            mean_fpr = np.linspace(0, 1, 100)
            fig, ax = plt.subplots()
            feature_importance = []
            for i, (train, test) in enumerate(cv.split(x, y)):
                clf = classifier.fit(x[train], y[train])

                viz = plot_roc_curve(classifier, x[test], y[test],
                                     name='ROC fold {}'.format(i),
                                     alpha=0.3, lw=1, ax=ax)
                interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
                interp_tpr[0] = 0.0
                tprs.append(interp_tpr)
                aucs.append(viz.roc_auc)
                print(i)

            ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
                    label='Chance', alpha=.8)

            mean_tpr = np.mean(tprs, axis=0)
            mean_tpr[-1] = 1.0
            mean_auc = auc(mean_fpr, mean_tpr)
            std_auc = np.std(aucs)
            ax.plot(mean_fpr, mean_tpr, color='b',
                    label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),lw=2, alpha=.8)

            std_tpr = np.std(tprs, axis=0)
            tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
            tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
            ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                            label=r'$\pm$ 1 std. dev.')

            ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05])
            ax.legend(loc="lower right")
            plt.savefig(output_folder + 'ROC_SVC_Adaboost_' + comparison + '_lr' +
                        str(lr) + '_n' + str(n) + '_k5.png', dpi=300)
            plt.clf()

    # using the best hyperparameters
    classifier = AdaBoostClassifier(base_estimator=svc, n_estimators=100,
                                    learning_rate=0.01, random_state=0)
    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)
    fig, ax = plt.subplots()
    feature_importance = []
    for i, (train, test) in enumerate(cv.split(x, y)):
        clf = classifier.fit(x[train], y[train])
        """For feature importance calculation"""
        results = permutation_importance(clf, x[test], y[test], random_state = 0)
        importance = results.importances_mean
        feature_importance.append(importance)

        """This is for drawing individual plots"""
        viz = plot_roc_curve(classifier, x[test], y[test],
                             name='Set {}'.format(i+1),
                             alpha=0.3, lw=1, ax=ax)
        interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
        interp_tpr[0] = 0.0
        tprs.append(interp_tpr)
        aucs.append(viz.roc_auc)
        print(i)
    
    """For feature importance calculation"""
    feature_importance = np.asarray(feature_importance).T
    mean_feature_importance = np.mean(feature_importance, axis=1)
    std_feature_importance = np.std(feature_importance, axis=1)
    df = pd.DataFrame({'Gene': gene_list, 'Importance_scores': list(feature_importance),
                       'Mean': list(mean_feature_importance), 'Std': list(std_feature_importance)})
    df.to_csv(output_folder + "FI_SVC_Adaboost_lr0.01_n100_k5_nr5.tsv", sep='\t', index=False)


    ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r', alpha=.8)

    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    ax.plot(mean_fpr, mean_tpr, color='b',
            label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),lw=2, alpha=.8)

    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                    label=r'$\pm$ 1 std. dev.')

    ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05])
    ax.legend(loc="lower right")

    plt.savefig(output_folder + 'FINALFIGURE_ROC_SVC_Adaboost_' + comparison + '_lr0.01_n100_k5.png', dpi=300)
