'''
Created on Jan 28, 2020
(Started pretty much re-writing on Apr 8, 2020)

@author: ywkim
'''

import sys
import csv
import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path
from random import shuffle
import matplotlib.pyplot as plt
from germline_analyses import get_matrix_as_df
from sklearn.linear_model import LinearRegression


def get_residual(a_, c_, x_, y_):
    # ax+by+c=0 (b=-1)
    # r = ax-y+c/sqrt(a^2+b^2)
    y_hat = (a_ * x_ + c_)
    r_ = y_ - y_hat

    return r_


def get_regression_plot(x_reg, y_reg, xx, yy, x_label, y_label, title, test_type):
    plt.axhline(0, color='0.75')
    plt.axvline(0, color='0.75')
    sns.regplot(x = np.asarray(xx), y=np.asarray(yy), fit_reg=False)
    ax = sns.regplot(x = np.asarray(x_reg), y=np.asarray(y_reg), scatter=False)
    ax.set(xlabel=x_label, ylabel=y_label)
    if y_label == 'sum of EA':
        plt.xlim(xmin=0)
        plt.ylim(ymin=0)
    else:
        pass
    plt.title(title)
    plt.savefig(output_folder + '/plots/' + title + test_type + '.png', dpi = 200)
    plt.clf()


def initial_lin_reg_one_fit(df_ADe2_sum_, df_ADe2_freq_, df_HCe4_sum_, df_HCe4_freq_, test_type):
    # First establishing the overall fit
    df_ADe2_sum_['Total'] = df_ADe2_sum_.sum(axis=1)
    df_HCe4_sum_['Total'] = df_HCe4_sum_.sum(axis=1)
    df_ADe2_freq_['Total'] = df_ADe2_freq_.sum(axis=1)
    df_HCe4_freq_['Total'] = df_HCe4_freq_.sum(axis=1)

    X_ADe2 = df_ADe2_freq_['Total'].values.tolist()
    X_HCe4 = df_HCe4_freq_['Total'].values.tolist()
    Y_ADe2 = df_ADe2_sum_['Total'].values.tolist()
    Y_HCe4 = df_HCe4_sum_['Total'].values.tolist()

    X_total = X_ADe2 + X_HCe4
    Y_total = Y_ADe2 + Y_HCe4

    X_mat = np.transpose(np.asarray[X_total])
    Y_mat = np.transpose(np.asarray[Y_total])
    lin_reg = LinearRegression()

    lin_reg.fit(X_mat, Y_mat)
    a = lin_reg.coef_[0][0]
    c = lin_reg.intercept_[0]

    residual_ADe2 = []
    residual_HCe4 = []

    for idx, gene in enumerate(total_gene_list):
        x_ADe2 = df_ADe2_freq_['Total'][idx]
        y_ADe2 = df_ADe2_sum_['Total'][idx]
        x_HCe4 = df_HCe4_freq_['Total'][idx]
        y_HCe4 = df_HCe4_sum_['Total'][idx]

        r_ADe2 = get_residual(a, c, x_ADe2, y_ADe2)
        r_HCe4 = get_residual(a, c, x_HCe4, y_HCe4)

        residual_ADe2.append(r_ADe2)
        residual_HCe4.append(r_HCe4)

    df_ADe2 = pd.DataFrame({'Gene':total_gene_list, 'var_count':X_ADe2, 'sum_EA':Y_ADe2, 'Residual':residual_ADe2})
    df_HCe4 = pd.DataFrame({'Gene':total_gene_list, 'var_count':X_HCe4, 'sum_EA':Y_HCe4, 'Residual':residual_HCe4})
    df_ADe2.to_csv(output_folder + 'initial_residuals_ADe2.tsv', sep='\t', index=False)
    df_HCe4.to_csv(output_folder + 'initial_residuals_HCe4.tsv', sep='\t', index=False)

    get_regression_plot(X_total, Y_total, X_ADe2, Y_ADe2, '# of variants', 'sum of EA', 'residuals_ADe2', test_type)
    get_regression_plot(X_total, Y_total, X_HCe4, Y_HCe4, '# of variants', 'sum of EA', 'residuals_HCe4', test_type)

    return residual_HCe4, residual_ADe2


def lin_reg_of_residuals(x_list, y_list, random_ideal_dict_, test_type):
    x_mat = np.transpose(np.asarray[x_list])
    y_mat = np.transpose(np.asarray[y_list])
    lin_reg = LinearRegression()

    lin_reg.fit(x_mat, y_mat)
    a = lin_reg.coef_[0][0]
    c = lin_reg.intercept_[0]

    r_of_r_list = []

    for idx_, gene_ in enumerate(total_gene_list):
        x = x_list(idx_)
        y = y_list(idx_)
        r_of_residuals = get_residual(a, c, x, y)

        if test_type == 'test':
            r_of_r_list.append(r_of_residuals)
        elif test_type == 'random':
            try:
                random_ideal_dict_[gene_].append(r_of_residuals)
            except KeyError:
                random_ideal_dict_[gene_] = []
                random_ideal_dict_[gene_].append(r_of_residuals)
        else:
            print('Error: you should never get this error')
            sys.exit()

    df_ = pd.DataFrame({'Gene': total_gene_list, 'r_HCe4': x_list, 'r_ADe2': y_list, 'iDEAL': r_of_r_list})
    df_.to_csv(output_folder + 'r_of_residuals_iDEAL.tsv', sep='\t', index=False)
    get_regression_plot(x_list, y_list, x_list, y_list, 'r_HCe4', 'r_ADe2', 'iDEAL', test_type)

    return r_of_r_list


if __name__ == '__main__':
    # input_folder = '/Users/ywkim/rosinante/ADSP/iDEAL_input_folder/' # This would change depending on
    input_folder = '/lab/rosinante/shared/ADSP/iDEAL_input_folder/'   # where you have Rosinante mounted on
    output_folder = str(Path().absolute()) + '/output/'  # For now the output files will be stored locally

    df_ADe2_sum = get_matrix_as_df(output_folder, 'sum_ea_matrix_ADe2.tsv', sep='\t', index_col=0)
    df_HCe4_sum = get_matrix_as_df(output_folder, 'sum_ea_matrix_HCe4.tsv', sep='\t', index_col=0)
    df_ADe2_freq = get_matrix_as_df(output_folder, 'frequency_matrix_ADe2.tsv', sep='\t', index_col=0)
    df_HCe4_freq = get_matrix_as_df(output_folder, 'frequency_matrix_HCe4.tsv', sep='\t', index_col=0)

    df_genes = pd.read_csv(output_folder + 'total_gene_list.csv', sep=',', index_col=None)
    total_gene_list = df_genes['Gene'].values.tolist()
    print(len(total_gene_list))

    random_ideal_dict = {}

    residuals_HCe4, residuals_ADe2 = initial_lin_reg_one_fit(df_ADe2_sum, df_ADe2_freq, df_HCe4_sum, df_HCe4_freq,
                                                             '_test')

    ideal_list = lin_reg_of_residuals(residuals_HCe4, residuals_ADe2, random_ideal_dict, '_test')

    """Starting randomization step"""
    pt_HCe4 = df_HCe4_sum.columns.tolist()
    pt_ADe2 = df_ADe2_sum.columns.tolist()
    total_pt = pt_ADe2 + pt_HCe4

    random_e4_total = []
    random_e2_total = []

    for p in range(1000):
        random_e4 = []
        random_e2 = []
        shuffle(total_pt)

        for i in range(len(total_pt)):
            pt = total_pt[i]  # pt = each patient

            if i < len(pt_HCe4):
                random_e4.append(pt)
            else:
                random_e2.append(pt)

        random_e4_total.append(random_e4)
        random_e2_total.append(random_e2)

    df = pd.DataFrame(random_e4_total)
    df.to_csv(output_folder + 'random_e4_pt_list.tsv', sep='\t', index=False)

    df = pd.DataFrame(random_e2_total)
    df.to_csv(output_folder + 'random_e2_pt_list.tsv', sep='\t', index=False)

    df_sum = pd.merge(df_HCe4_sum, df_ADe2_sum, left_index=True, right_index=True)  # default = inner join
    df_freq = pd.merge(df_HCe4_freq, df_ADe2_freq, left_index=True, right_index=True)

    for pt_idx in range(1000):
        random_e4 = random_e4_total[pt_idx]
        random_e2 = random_e2_total[pt_idx]

        df_random_e4_sum = df_sum[random_e4]
        df_random_e2_sum = df_sum[random_e2]
        df_random_e4_freq = df_freq[random_e4]
        df_random_e2_freq = df_freq[random_e2]

        residuals_random_e4, residuals_random_e2 = initial_lin_reg_one_fit(df_random_e2_sum,
                                                                           df_random_e2_freq,
                                                                           df_random_e4_sum,
                                                                           df_random_e4_freq, '_random')

        lin_reg_of_residuals(residuals_random_e4, residuals_random_e2, random_ideal_dict, '_random')

        # TODO save ideal_random_dict

    # Get z-scores

    outputFile = output_folder + 'iDEAL_with_z-scores.tsv'
    with open(outputFile, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Gene', 'iDEAL', 'mean', 'std', 'z-score', 'iDEAL_random'])
        for idx, gene in total_gene_list:
            random_ideal_list = random_ideal_dict[gene]
            mean = np.mean(random_ideal_list)
            std = np.std(random_ideal_list)
            xr = ideal_list[idx]
            z = (xr - mean) / std
            # new_r
            # for d in random_ideal_list:
            #    new_r.append(float('{0:.2f}'.format(d)))
            # info = ([gene, float('{0:.2f}'.format(x)), float('{0:.2f}'.format(mean)),
            #          float('{0:.2f}'.format(std)), float('{0:.2f}'.format(z)), new_r])

            info = ([gene, xr, mean, std, z, random_ideal_list])
            writer.writerow(info)
