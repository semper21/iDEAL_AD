'''
Created on Mar 3, 2020

@author: ywkim
'''

import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from collections import Counter

from germline_analyses import matrix_out, output_dict

def get_quality(input_file):
    dict_ = {}
    for line in open(input_file):
        if line[0] == '#':
            continue
        cols = line.strip().split('\t')
        chromosome = cols[0]
        pos = cols[1]
        ref = cols[2]
        alt = cols[3]

        filter_ = cols[4]

        s = '-'
        mut = (chromosome, pos, ref, alt)
        substitution = s.join(mut)

        dict_[substitution] = filter_

    # Saving the dictionary to a local file since it will be used in other codes as well
    pd.DataFrame.from_dict(dict_, orient='index').to_csv(input_folder + 'ADSP_quality_file.csv')

    return dict_


def get_phenotype(pheno_file):
    AD2, AD3, AD4, HC2, HC3, HC4 = [], [], [], [], [], []
    for line in open(pheno_file):
        if line[0] == '#' or line[0] == 'd':
            continue
        cols = line.strip().split('\t')

        sub_id = cols[1]
        apoe = str(cols[7])
        race = str(cols[10])
        eth = str(cols[11])
        state = cols[13]

        if race == '5' and eth == '0':  # Caucasians only
            if str(state) == '1':  # AD
                if apoe == '22' or apoe == '23':  # APOE2
                    AD2.append(sub_id)
                elif apoe == '33':  # APOE3
                    AD3.append(sub_id)
                elif apoe == '44' or apoe == '34':  # APOE4
                    AD4.append(sub_id)  # healthy with risk
                else:
                    pass

            elif str(state) == '0':  # HC
                if apoe == '22' or apoe == '23':  # APOE2
                    HC2.append(sub_id)
                elif apoe == '33':  # APOE3
                    HC3.append(sub_id)
                elif apoe == '44' or apoe == '34':  # APOE4
                    HC4.append(sub_id)
                else:
                    pass
            else:
                pass

    return AD2, AD3, AD4, HC2, HC3, HC4


def get_phenotype_from_excel(phenotype_file):
    df_= pd.read_excel(phenotype_file)
    HC2 = df_.loc[(df_['AD'] == 0) & (df_['APOE'].isin([22, 23]))]['ID'].values.tolist()
    HC3 = df_.loc[(df_['AD'] == 0) & (df_['APOE'] == 33)]['ID'].values.tolist()
    HC4 = df_.loc[(df_['AD'] == 0) & (df_['APOE'].isin([44, 34]))]['ID'].values.tolist()
    AD2 = df_.loc[(df_['AD'] == 1) & (df_['APOE'].isin([22, 23]))]['ID'].values.tolist()
    AD3 = df_.loc[(df_['AD'] == 1) & (df_['APOE'] == 33)]['ID'].values.tolist()
    AD4 = df_.loc[(df_['AD'] == 1) & (df_['APOE'].isin([44, 34]))]['ID'].values.tolist()

    return AD2, AD3, AD4, HC2, HC3, HC4


def string_to_ea(current_list):
    temp_list = []
    for ea in current_list:
        try:
            temp_list.append(float(ea))
        except ValueError:  # When there is a string
            if ea in ['STOP', 'no_STOP', 'STOP-loss', 'START_loss']:
                temp_list.append(1.0)
            else:
                continue
    return temp_list


def fill_matrix(trauma_file, ptidx, cohort):
    """

    :param trauma_file: each trauma file contains germline variant info for each patient
    :param ptidx: patient indices
    :return: fills in the matrices
    """
    ea_dict = {}
    variant_counter = Counter()
    for line in open(trauma_file):
        if line[0] == '#':
            continue
        cols = line.strip().split('\t')

        chro = cols[0]
        pos = cols[1]
        ref = cols[2]
        alt = cols[3]
        gene = cols[4].split(';')[0]
        sub = cols[5]
        action = cols[6]

        if gene == '':
            continue
        if cohort == 'ADSP_discovery':
            # Only needed for APSP_discovery
            s = '-'
            mut = (chro, pos, ref, alt)
            substitution = s.join(mut)

            try:
                qual = dict_qc[substitution]
            except KeyError:
                qual = 'non-PASS'

            if qual == 'non-PASS':
                continue
            elif qual == 'PASS':
                pass
            else:
                print('ERROR')

        if sub == '-':  # no indels
            continue
        if gene == 'APOE' and sub == 'C130R':
            continue

        # get number of all variants (no maf cutoff, all var(except indels))
        variant_counter[gene] += 1

        # get sum EA
        if gene not in ea_dict:
            ea_dict[gene] = []
        if action in ['silent', '-', 'no_action', 'no_trace', 'no_gene'] or gene == '':
            continue
        elif action in ['STOP', 'no_STOP', 'STOP-loss', 'START_loss']:
            ea_dict[gene].append(1.)
        else:
            try:
                spl_ea = action.strip().split(';')
                new_spl_ea = string_to_ea(spl_ea)
                if not new_spl_ea:
                    continue
                else:
                    average = np.mean(new_spl_ea)
                    ea_dict[gene].append(average / 100.)
            except AttributeError:
                ea_dict[gene].append(float(action) / 100.)

    for gene_idx, gene in enumerate(total_gene_list):
        try:
            ea_list = ea_dict[gene]
            sum_ea = np.sum(ea_list)
        except KeyError:
            sum_ea = 0

        try:
            freq = variant_counter[gene]
        except KeyError:
            freq = 0

        matrix_sum[gene_idx][ptidx] = sum_ea
        matrix_freq[gene_idx][ptidx] = freq


if __name__ == '__main__':
    cohort_name = 'ADSP_discovery'
    # cohort_name = 'ADSP_extension'
    input_folder = '/lab/rosinante/shared/ADSP/iDEAL_input_folder/' + cohort_name + '/'  # where Rosinante is mounted on

    output_folder = str(Path().absolute()) + '/output_' + cohort_name + '/' # for now the output files
                                                                            # will be stored locally
    # FOR ADSP_DISCOVERY
    if cohort_name == 'ADSP_discovery':
        quality_file = input_folder + 'snvquality_detailed_jamie.csv'  # generated using the raw data from dbgap
        dict_qc = get_quality(quality_file)

        phenotype_file = input_folder + 'phs000572.v7.pht005179.v1.p4.c1.CaseControlEnrichedPhenotypesWES_y1.HMB-IRB.txt'
        ADe2, ADe3, ADe4, HCe2, HCe3, HCe4 = get_phenotype(phenotype_file)

        trauma_folder = input_folder + 'all_short_name/'

    # FOR ADSP_EXTENSION
    elif cohort_name == 'ADSP_extension':
        phenotype_file = input_folder + 'ADSP_extension_phenotypes.xlsx'
        ADe2, ADe3, ADe4, HCe2, HCe3, HCe4 = get_phenotype_from_excel(phenotype_file)

        trauma_folder = input_folder + 'Actions/'

    else:
        print('input folder invalid')
        sys.exit()

    df_genes = pd.read_csv(output_folder + 'total_gene_list.csv', sep=',', index_col=None)
    total_gene_list = df_genes['Gene'].values.tolist()
    print(len(total_gene_list))

    pt_lists = [ADe2, ADe3, ADe4, HCe2, HCe3, HCe4]
    for idx, group in enumerate(['ADe2', 'ADe3', 'ADe4', 'HCe2', 'HCe3', 'HCe4']):
        pt_list = pt_lists[idx]
        matrix_freq = np.zeros((len(total_gene_list), len(pt_list)))
        matrix_sum = np.zeros((len(total_gene_list), len(pt_list)))
        counter = 0
        for pt_idx, pt_id in enumerate(pt_list):
            if cohort_name == 'ADSP_discovery':
                pt_file = os.path.join(trauma_folder, pt_id)
            elif cohort_name == 'ADSP_extension':
                pt_file = os.path.join(trauma_folder, pt_id + '.trauma')
            else:
                print('input folder invalid')
                sys.exit()
            fill_matrix(pt_file, pt_idx, cohort_name)
            counter += 1

            print(group, counter)

        matrix_out(matrix_freq, pt_list, total_gene_list, output_folder, 'frequency_matrix_' + group)
        matrix_out(matrix_sum, pt_list, total_gene_list, output_folder, 'sum_ea_matrix_' + group)
        output_dict({k: v for v, k in enumerate(pt_list)}, output_folder, group + '_ptidx.tsv', sep='\t')