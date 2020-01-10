"""
Created on Aug 23, 2019

@author: ywkim

10/24/29
-Output file is more cleaned up
-Action: average (if more than one per variant)
-CI format: x.xx
"""

import csv
import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels.api as sm

from IPython import embed


def calculate_or(a, b, c, d):
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


def calculate_or_with_ci(a, b, c, d):
    table = np.zeros((2, 2))
    # a: affected with mutation
    table[0][0] = a
    # b: healthy with mutation
    table[0][1] = b
    # c: affected w/o mutation
    table[1][0] = c
    # d: healthy w/o mutation
    table[1][1] = d

    or_table = sm.stats.Table2x2(table)

    odds_ratio = or_table.oddsratio
    p_value = or_table.oddsratio_pvalue()
    confidence_interval = list(or_table.oddsratio_confint())

    '''used for cleaned-up version'''
    confidence_interval = [float(f'{x:.2f}') for x in confidence_interval]

    return odds_ratio, p_value, confidence_interval


if __name__ == '__main__':
    dirc = '/Users/ywkim/Desktop/Projects/GermlineProject/ADSP/RVEA/new2018/' \
           'RVEA_BaylorPass_nonHisWhite_2v4_h5py_STARTLOSS100/'

    # genes = ['LRRC17', 'TRAF3IP', 'UGT3A1', 'ARHGAP44', 'PRKAG3', 'THG1L', 'GPR37', 'ATP6V0E2', 'C1orf185']

    '''
    fly_hits_file = dirc + 'hit_genes.txt'
    df_hits = pd.read_csv(fly_hits_file)
    genes = df_hits['Gene'].values.tolist()
    '''

    # all 216 genes
    ideal_file = dirc + 'iDEAL_genelist.txt'
    df_genes = pd.read_csv(ideal_file)
    genes = df_genes['Gene'].values.tolist()

    caseFile = dirc + '216_updated/216_mutation_count_APOE3-AD.csv'
    controlFile = dirc + '216_updated/216_mutation_count_APOE3-HC.csv'

    dfCase = pd.read_csv(caseFile, sep=',')
    #dfCase = dfCase[['Gene', 'Sub', 'Action', 'Count', 'Het:Hom']]
    dfControl = pd.read_csv(controlFile, sep=',')
    #dfControl = dfControl[['Gene', 'Sub', 'Action', 'Count', 'Het:Hom']]

    df = dfCase.merge(dfControl, on=['Gene', 'Sub', 'Action'], how='outer')
    df = df.fillna(0)
    df.columns = ['Gene', 'Sub', 'Action', '# APOE3-AD', 'Het APOE3-AD', 'Hom APOE3-AD', 'Het:Hom APOE3-AD',
                  '# APOE3-HC', 'Het APOE3-HC', 'Hom APOE3-HC', 'Het:Hom APOE3-HC']

    #df.to_csv(dirc + '216_updated/216_all_mutation_count_APOE3_withEA.csv', sep=',', index=False)

    hc_all = 1657   # only 3/3 HC
    ad_all = 1346   # only 3/3 AD

    outputFile = dirc + '216_updated/216genes_mutation_count_APOE3_withEA_OR_CI_clean2.csv'
    #outputFile = dirc + '216_updated/69hits_all_mutation_count_APOE3_withEA_OR_CI_clean.csv'
    with open(outputFile, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(['Gene', 'Sub', 'Action', 'Het:Hom Case', 'Het:Hom Control', 'Het_OR', 'Het_CI', 'Het_p-value',
                         'Hom_OR', 'Hom_CI', 'Hom_p-value', '# Case', '# Control', 'OR', 'CI', 'p-value'])
        for gene in genes:
            sub_info = df[df['Gene'] == gene].values
            for row in sub_info:
                sub = row[1]
                action = row[2]
                ad = row[3]
                het_ad = row[4]
                hom_ad = row[5]
                ratio_ad = row[6]
                hc = row[7]
                het_hc = row[8]
                hom_hc = row[9]
                ratio_hc = row[10]

                diff = int(ad) - int(hc)

                if action in ['silent', 'no_action', 'no_gene', 'no_trace']:
                    continue

                else:
                    if ad == 0 or hc == 0:
                        odds_ratio, p_value, ci = calculate_or_with_ci(ad + 1, hc + 1, ad_all - ad + 1, hc_all - hc + 1)
                    else:
                        odds_ratio, p_value, ci = calculate_or_with_ci(ad, hc, ad_all - ad, hc_all - hc)

                    if het_ad == 0 and het_hc == 0:
                        het_or = 'NA'
                        het_ci = 'NA'
                        het_p = 'NA'
                    elif het_ad == 0 or het_hc == 0:
                        het_or, het_p, het_ci = calculate_or_with_ci(het_ad + 1, het_hc + 1, ad_all - het_ad + 1, hc_all - het_hc + 1)
                    else:
                        het_or, het_p, het_ci = calculate_or_with_ci(het_ad, het_hc, ad_all - het_ad, hc_all - het_hc)

                    if hom_ad == 0 and hom_hc == 0:
                        hom_or = 'NA'
                        hom_ci = 'NA'
                        hom_p = 'NA'
                    elif hom_ad == 0 or hom_hc == 0:
                        hom_or, hom_p, hom_ci = calculate_or_with_ci(hom_ad + 1, hom_hc + 1, ad_all - hom_ad + 1, hc_all - hom_hc + 1)
                    else:
                        hom_or, hom_p, hom_ci = calculate_or_with_ci(hom_ad, hom_hc, ad_all - hom_ad, hc_all - hom_hc)

                    if ratio_ad == 0:
                        ratio_ad = '0:0'
                    if ratio_hc == 0:
                        ratio_hc = '0:0'

                    '''used for cleaned-up version'''
                    actions = action.split(';')
                    if len(actions) == 1:
                        mean_action = action
                    else:
                        temp_list = []
                        for x in actions:
                            try:
                                temp_list.append(float(x))
                            except ValueError: # this is for when Action contains 'no_gene'
                                pass
                        mean_action = np.mean(temp_list)

                    info = [gene, sub, mean_action, ratio_ad, ratio_hc, het_or, het_ci, het_p, hom_or, hom_ci, hom_p,
                            ad, hc, odds_ratio, ci, p_value]
                    writer.writerow(info)
