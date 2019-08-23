'''
Created on 

@author: ywkim
'''
import csv
import numpy as np
import pandas as pd
import scipy.stats as stats


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


if __name__ == '__main__':
    dirc = '/Users/ywkim/Desktop/Projects/GermlineProject/ADSP/RVEA/new2018/' \
           'RVEA_BaylorPass_nonHisWhite_2v4_h5py_STARTLOSS100/'

    fly_hits_file = dirc + 'hit_genes.txt'
    df_hits = pd.read_csv(fly_hits_file)
    genes = df_hits['Gene'].values.tolist()

    case_file = dirc + '216_updated/216_mutation_count_APOE2-AD.csv'
    control_file = dirc + '216_updated/216_mutation_count_APOE4-HC.csv'

    df_case = pd.read_csv(case_file, sep=',')
    df_case = df_case[['Gene', 'Sub', 'Action', 'Count', 'Het:Hom']]
    df_control = pd.read_csv(control_file, sep=',')
    df_control = df_control[['Gene', 'Sub', 'Action', 'Count', 'Het:Hom']]

    df = df_case.merge(df_control, on=['Gene', 'Sub', 'Action'], how='outer')
    df = df.fillna(0)
    df.columns = ['Gene', 'Sub', 'Action', '# APOE2-AD', 'Het:Hom APOE2-AD', '# APOE4-HC', 'Het:Hom APOE4-HC']
    df.to_csv(dirc + '216_all_mutation_count_APOE3_withEA.csv', sep=',', index=False)
    df.to_csv(dirc + '216_updated/216_all_mutation_count_APOE2v4_withEA.csv', sep=',', index=False)

    hc_all = 301
    ad_all = 179

    # only 69 genes
    outputFile = dirc + '216_updated/69hits_all_mutation_count_APOE2v4_withEA_OR.csv'
    with open(outputFile, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(['Gene', 'Sub', 'Action', 'Het:Hom Case', 'Het:Hom Control', '# Case', '# Control', 'OR', 'p-value'])

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
                        oddsratio, pvalue = calculate_or(ad + 1, hc + 1, ad_all - ad + 1, hc_all - hc + 1)
                    else:
                        oddsratio, pvalue = calculate_or(ad, hc, ad_all - ad, hc_all - hc)

                    info = [gene, sub, action, ratio_ad, ratio_hc, ad, hc, oddsratio, pvalue]
                    writer.writerow(info)
