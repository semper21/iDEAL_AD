"""


@author: ywkim
"""

import csv
import numpy as np
import pandas as pd
from sys import argv
import seaborn as sns
from scipy import stats
import matplotlib.pyplot as plt


def getCorr(listA, listB, name):
    sns.set_style("whitegrid", {'axes.grid': False})
    sns.color_palette("cubehelix", 2)
    plt.axhline(0, color='0.75')
    plt.axvline(0, color='0.75')

    slope, intercept, r_value, p_valueS, std_err = stats.linregress(listA, listB)

    sns.regplot(x=np.asarray(listA), y=np.asarray(listB), fit_reg=False, scatter_kws={"color": "0.1", "s": 16})
    ax = sns.regplot(x=np.asarray(listA), y=np.asarray(listB), scatter=False, color="darkred")
    ax.set(xlabel='Ti/Tv', ylabel='lambda')
    plt.xlim(xmin=1.)
    plt.ylim(ymin=0.01)
    plt.savefig(dirc + 'LambdaVsTiTv_' + name + '.png', dpi=300)
    plt.clf()

    rsquared = r_value ** 2
    cc, p_valueP = stats.pearsonr(listA, listB)

    outputFile = (dirc + 'LambdaVsTiTv_' + name + '.csv')

    with open(outputFile, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['r:', r_value])
        writer.writerow(['p-value:', p_valueS])
        writer.writerow(['r-squared:', rsquared])
        writer.writerow(['pearson corr coeff:', cc])
        writer.writerow(['p-value:', p_valueP])

if __name__ == '__main__':
    dirc = argv[1]      # 'Quality' folder

    lambda_total = []
    titv_action_total = []
    titv_all_total = []

    outputFile = dirc + 'percent_false_positives.csv'
    with open(outputFile, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        filename = dirc + 'after_QCfilter_jamie.exome_qual'
        df = pd.read_csv(filename, sep='\t', index_col=False)
        lambda_list = df['Lambda'].values.tolist()
        titv_action_list = df['TiTv_act'].values.tolist()
        titv_all_list = df['TiTv_all'].values.tolist()
        lambda_total += lambda_list
        titv_action_total += titv_action_list
        titv_all_total += titv_all_list

        noise_list = []
        for val in lambda_list:
            l = float(val) / 0.03752
            noise = (np.log(l) / (-0.01261))
            if noise < 0:
                noise = 0
            elif noise > 100:
                noise = 100
            else:
                pass
            noise_list.append(noise)
        # print len(newlistL)
        info = ['ADSP'] + noise_list
        writer.writerow(info)

    # print listAllL
    # print listAllT
    getCorr(titv_action_total, lambda_total, 'missense')
    getCorr(titv_all_total, lambda_total, 'allcoding')

    with open(dirc + 'averageValues.csv', 'w') as output:
        writer = csv.writer(output, lineterminator='\n')
        writer.writerow(['average lambda = ', np.mean(lambda_total), '+-', np.std(lambda_total)])
        writer.writerow(['median lambda = ', np.median(lambda_total), '+-', np.std(lambda_total)])
        writer.writerow(['average TiTv = ', np.mean(titv_all_total), '+-', np.std(titv_all_total)])
        writer.writerow(['median TiTv = ', np.median(titv_all_total), '+-', np.std(titv_all_total)])
    print('average lambda = ', np.mean(lambda_total), '+-', np.std(lambda_total))
    print('average TiTv = ', np.mean(titv_all_total), '+-', np.std(titv_all_total))
