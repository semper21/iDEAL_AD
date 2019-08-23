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
           'RVEA_BaylorPass_nonHisWhite_2v4_h5py_STARTLOSS100/216genes_variants/'

    caseFile = dirc + '216_mutation_count_APOE3-AD.txt'
    controlFile = dirc + '216_mutation_count_APOE3-HC.txt'

    dfCase = pd.read_csv(caseFile, sep=',')
    dfControl = pd.read_csv(controlFile, sep=',')

    df = dfCase.merge(dfControl, on=['Gene', 'Sub'], how='outer')
    df = df.fillna(0)

    df.to_csv(dirc + 'hippoGenes_all_mutation_count_withEA.txt', sep='\t', index=False)

    hc_all = 1657   # only 3/3 HC
    ad_all = 1346   # only 3/3 AD

    outputFile = dirc + '216_all_mutation_count_APOE3_withEA_OR.txt'
    with open(outputFile, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Gene', 'Sub', 'Action', 'Case', 'Control', 'OR', 'p-value'])

        for gene in genes:
            sub_info = df[df['Gene'] == gene].values
            for row in sub_info:

                mut = row[1]
                ad = row[2]
                hc = row[3]

                subs = mut.split('(')
                sub = subs[0]
                action = ''.join(list(subs[1])[:-1])

                sub_one = sub.split(';')[0]
                action_one = action.split(';')[0]

                diff = int(ad) - int(hc)

                if action in ['silent', 'no_action']:
                    continue

                else:
                    # TODO : change the code to just calculate OR as is
                    '''
                    if diff > 0: #pathogenic
                        if ad == 0 or hc == 0:
                            oddsratio, pvalue = calcOR(ad+1, hc+1, ad_all-ad+1, hc_all-hc+1)
                        else:
                            oddsratio, pvalue = calcOR(ad, hc, ad_all-ad, hc_all-hc)


                    elif diff < 0: #protective
                        if ad == 0 or hc == 0:
                            oddsratio, pvalue = calcOR(hc+1, ad+1, hc_all-hc+1, ad_all-ad+1)
                        else:
                            oddsratio, pvalue = calcOR(hc, ad, hc_all-hc, ad_all-ad)
                    else:
                        continue

                    '''

                    if ad == 0 or hc == 0:
                        oddsratio, pvalue = calcOR(ad + 1, hc + 1, ad_all - ad + 1, hc_all - hc + 1)
                    else:
                        oddsratio, pvalue = calcOR(ad, hc, ad_all - ad, hc_all - hc)
                    # '''

                    info = [gene, sub, action, ad, hc, oddsratio, pvalue]
                    writer.writerow(info)
