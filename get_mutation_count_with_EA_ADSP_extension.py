"""
Created on 

@author: ywkim
"""

import os
import csv
import sys
import numpy as np
import pandas as pd
import seaborn as sns
from collections import Counter
import matplotlib.pyplot as plt
from pathlib import Path
from IPython import embed


def get_phenotype_from_excel(phenotype_file):
    df_= pd.read_excel(phenotype_file)
    HC2 = df_.loc[(df_['AD'] == 0) & (df_['APOE'].isin([22, 23]))]['ID'].values.tolist()
    HC3 = df_.loc[(df_['AD'] == 0) & (df_['APOE'] == 33)]['ID'].values.tolist()
    HC4 = df_.loc[(df_['AD'] == 0) & (df_['APOE'].isin([44, 34]))]['ID'].values.tolist()
    AD2 = df_.loc[(df_['AD'] == 1) & (df_['APOE'].isin([22, 23]))]['ID'].values.tolist()
    AD3 = df_.loc[(df_['AD'] == 1) & (df_['APOE'] == 33)]['ID'].values.tolist()
    AD4 = df_.loc[(df_['AD'] == 1) & (df_['APOE'].isin([44, 34]))]['ID'].values.tolist()

    return AD2, AD3, AD4, HC2, HC3, HC4


def delNoStr(currList):
    tempList = []
    for x in currList:
        try:
            tempList.append(float(x))
        except ValueError:  # When There is a string
            if x in ['STOP', 'no_STOP', 'STOP-loss', 'START_loss']:
                tempList.append(1.)
            else:
                continue
    return tempList


def getInfo(ptFile, EA_dict, sub_counter, het_counter, homo_counter):
    for line in open(ptFile):
        if line[0] == '#':
            continue
        cols = line.strip().split('\t')
        chr = cols[0]
        pos = cols[1]
        ref = cols[2]
        alt = cols[3]
        gene = cols[4].split(';')[0]
        sub = cols[5]
        action = cols[6]
        zygosity = cols[7]

        if gene == '':
            continue

        s = '-'
        mut = (chr, pos, ref, alt)
        substitutionG = s.join(mut)

        if sub == '-':  # no indels
            continue
        if gene == 'APOE' and sub == 'C130R':
            continue

        # get sumEA
        if gene in gene_list:
            seq = sub + '(' + str(action) + ')'

            if gene not in EA_dict:
                EA_dict[gene] = []

            if action in ['silent', '-', 'no_action', 'no_trace', 'no_gene']:
                pass
            elif action in ['STOP', 'no_STOP', 'STOP-loss', 'START_loss']:
                EA_dict[gene].append(1.)
            else:
                try:
                    spl_EA = action.strip().split(';')
                    new_spl_EA = delNoStr(spl_EA)
                    if new_spl_EA == []:
                        pass
                    else:
                        average = np.mean(new_spl_EA)
                        EA_dict[gene].append(average / 100.)

                except AttributeError:
                    EA_dict[gene].append(float(action) / 100.)

            try:
                sub_counter[gene][seq] += 1

            except KeyError:
                sub_counter[gene] = Counter()
                sub_counter[gene][seq] += 1

            if str(zygosity) == '1':
                try:
                    het_counter[gene][seq] += 1
                except KeyError:
                    het_counter[gene] = Counter()
                    het_counter[gene][seq] += 1

            elif str(zygosity) == '2':
                try:
                    homo_counter[gene][seq] += 1
                except KeyError:
                    homo_counter[gene] = Counter()
                    homo_counter[gene][seq] += 1
            else:
                print('neither hetero-/homozygous')
                sys.exit()


def output_dict(sub_dict, het_dict, homo_dict, name):
    output_file = output_folder + '216_mutation_count_' + name + '.csv'
    with open(output_file, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(['Gene', 'Sub', 'Action', 'Count', 'Het_count', 'Hom_count', 'Het:Hom'])
        for gene in gene_list:
            try:
                sub_list = sub_dict[gene].keys()
                for sub_info in sub_list:
                    subs = sub_info.split('(')
                    sub = subs[0]
                    action = ''.join(list(subs[1])[:-1])

                    count = sub_dict[gene][sub_info]
                    try:
                        het_count = het_dict[gene][sub_info]
                    except KeyError: # if the mutation does not exist in heterozygous
                        het_count = 0

                    try:
                        homo_count = homo_dict[gene][sub_info]
                    except KeyError: # if the mutation does not exist in homozygous
                        homo_count = 0

                    ratio = str(het_count) + ':' + str(homo_count)
                    info = ([gene, sub, action, count, het_count, homo_count, ratio])
                    writer.writerow(info)

            except KeyError: # if the gene/mutation does not exist
                continue


if __name__ == '__main__':
    cohort_name = 'ADSP_extension'
    input_folder = '/lab/rosinante/shared/ADSP/iDEAL_input_folder/' + cohort_name + '/'  # where Rosinante is mounted on

    output_folder = str(Path().absolute()) + '/output_' + cohort_name + '/'  # for now the output files will be stored locally

    phenotype_file = input_folder + 'ADSP_extension_phenotypes.xlsx'
    apoe2_ad, apoe3_ad, ADe4, HCe2, apoe3_hc, apoe4_hc = get_phenotype_from_excel(phenotype_file)

    trauma_folder = input_folder + 'Actions/'

    ea_apoe2_ad = {}
    ea_apoe4_hc = {}
    ea_apoe3_ad = {}
    ea_apoe3_hc = {}

    sub_apoe2_ad = {}
    sub_apoe4_hc = {}
    sub_apoe3_ad = {}
    sub_apoe3_hc = {}

    het_count_apoe2_ad = {}
    het_count_apoe4_hc = {}
    het_count_apoe3_ad = {}
    het_count_apoe3_hc = {}

    homo_count_apoe2_ad = {}
    homo_count_apoe4_hc = {}
    homo_count_apoe3_ad = {}
    homo_count_apoe3_hc = {}

    gene_file = output_folder + 'iDEAL_og_genelist.txt'
    dfGene = pd.read_csv(gene_file, index_col=None)

    gene_list = dfGene['Gene'].values.tolist()

    print('-----Getting info-----')

    errorlog = []


    pt_count = 0
    for idx, pt_id in enumerate(apoe4_hc):
        pt_file = (os.path.join(trauma_folder, pt_id + '.trauma'))
        getInfo(pt_file, ea_apoe4_hc, sub_apoe4_hc, het_count_apoe4_hc, homo_count_apoe4_hc)
        pt_count += 1

    for idx, pt_id in enumerate(apoe3_hc):
        pt_file = (os.path.join(trauma_folder, pt_id + '.trauma'))
        getInfo(pt_file, ea_apoe3_hc, sub_apoe3_hc, het_count_apoe3_hc, homo_count_apoe3_hc)
        pt_count += 1

    print(pt_count)

    for idx, pt_id in enumerate(apoe2_ad):
        pt_file = (os.path.join(trauma_folder, pt_id + '.trauma'))
        getInfo(pt_file, ea_apoe2_ad, sub_apoe2_ad, het_count_apoe2_ad, homo_count_apoe2_ad)
        pt_count += 1

    for idx, pt_id in enumerate(apoe3_ad):
        pt_file = (os.path.join(trauma_folder, pt_id + '.trauma'))
        getInfo(pt_file, ea_apoe3_ad, sub_apoe3_ad, het_count_apoe3_ad, homo_count_apoe3_ad)

    print(pt_count)

    output_dict(sub_apoe2_ad, het_count_apoe2_ad, homo_count_apoe2_ad, 'APOE2-AD')
    output_dict(sub_apoe4_hc, het_count_apoe4_hc, homo_count_apoe4_hc, 'APOE4-HC')
    output_dict(sub_apoe3_ad, het_count_apoe3_ad, homo_count_apoe3_ad, 'APOE3-AD')
    output_dict(sub_apoe3_hc, het_count_apoe3_hc, homo_count_apoe3_hc, 'APOE3-HC')
