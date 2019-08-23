"""
Created on 

@author: ywkim
"""

import csv
import numpy as np
import pandas as pd
import scipy.stats as stats

from IPython import embed


def calcOR(a, b, c, d):
    table = np.zeros((2, 2))
    # a: affected with mutation
    table[0][0] = a
    # b: healthy with mutation
    table[0][1] = b
    # c: affected w/o mutation
    table[1][0] = c
    # d: healthy w/o mutation
    table[1][1] = d
    try:
        oddsratio, pvalue = stats.fisher_exact(table)
    except ValueError:
        embed()
    return oddsratio, pvalue


if __name__ == '__main__':
    genes = ['LRRC17', 'TRAF3IP', 'UGT3A1', 'ARHGAP44', 'PRKAG3', 'THG1L', 'GPR37', 'ATP6V0E2', 'C1orf185']

    dirc = '/Users/ywkim/Desktop/Projects/GermlineProject/ADSP/RVEA/new2018/' \
           'RVEA_BaylorPass_nonHisWhite_2v4_h5py_STARTLOSS100/216_updated/'

    caseFile = dirc + '216_mutation_count_APOE3-AD.csv'
    controlFile = dirc + '216_mutation_count_APOE3-HC.csv'

    dfCase = pd.read_csv(caseFile, sep=',')
    dfCase = dfCase[['Gene', 'Sub', 'Action', 'Count', 'Het:Hom']]
    dfControl = pd.read_csv(controlFile, sep=',')
    dfControl = dfControl[['Gene', 'Sub', 'Action', 'Count', 'Het:Hom']]

    df = dfCase.merge(dfControl, on=['Gene', 'Sub', 'Action'], how='outer')
    df = df.fillna(0)
    df.columns = ['Gene', 'Sub', 'Action', '# APOE3-AD', '# APOE3-HC', 'Het:Hom APOE3-AD', 'Het:Hom APOE3-HC']
    df.to_csv(dirc + '216_all_mutation_count_APOE3_withEA.csv', sep=',', index=False)

    hc_all = 1657   # only 3/3 HC
    ad_all = 1346   # only 3/3 AD

    outputFile = dirc + '216_all_mutation_count_APOE3_withEA_OR.csv'
    with open(outputFile, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(['Gene', 'Sub', 'Action', 'Het:Hom Case', 'Het:Hom Contrl', '# Case', '# Control', 'OR', 'p-value'])

        for gene in genes:
            sub_info = df[df['Gene'] == gene].values
            for row in sub_info:
                sub = row[1]
                action = row[2]
                ad = row[3]
                ratio_ad = row[4]
                hc = row[5]
                ratio_hc = row[6]
                diff = int(ad) - int(hc)

                if action in ['silent', 'no_action']:
                    continue

                else:
                    if ad == 0 or hc == 0:
                        oddsratio, pvalue = calcOR(ad + 1, hc + 1, ad_all - ad + 1, hc_all - hc + 1)
                    else:
                        oddsratio, pvalue = calcOR(ad, hc, ad_all - ad, hc_all - hc)

                    info = [gene, sub, action, ratio_ad, ratio_hc, ad, hc, oddsratio, pvalue]
                    writer.writerow(info)
