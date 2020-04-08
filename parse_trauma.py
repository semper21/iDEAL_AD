'''
Created on Mar 3, 2020

@author: ywkim
'''

import os
import numpy as np
import pandas as pd
from collections import Counter

from germline_analyses import delNoStr, output_dict, matrix_out

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
    for line in open(pheno_file):
        if line[0] == '#' or line[0] == 'd':
            continue
        cols = line.strip().split('\t')

        sub_id = cols[1]
        apoe = str(cols[7])
        race = str(cols[10])
        eth = str(cols[11])
        state = cols[13]

        AD2, AD3, AD4, HC2, HC3, HC4 = [], [], [], [], [], []

        if race == '5' and eth == '0':  # Caucasians only
            if str(state) == '1':                   # AD
                if apoe == '22' or apoe == '23':    # APOE2
                    AD2.append(sub_id)
                elif apoe == '33':                  # APOE3
                    AD3.append(sub_id)
                elif apoe == '44' or apoe == '34':  # APOE4
                    AD4.append(sub_id)    # healthy with risk
                else:
                    pass

            elif str(state) == '0':                 # HC
                if apoe == '22' or apoe == '23':    # APOE2
                    HC2.append(sub_id)
                elif apoe == '33':                  # APOE3
                    HC3.append(sub_id)
                elif apoe == '44' or apoe == '34':  # APOE4
                    HC4.append(sub_id)
                else:
                    pass
            else:
                pass

    return AD2, AD3, AD4, HC2, HC3, HC4


def fill_matrix(trauma_file, pt_idx):
    """ Fill matrices with average EA per gene per patient

    :param trauma_file: trauma file for an individual patient
    :param ptidx: index of each patient
    :return: Fills in 1) variant count matrix, 2) sum EA matrix
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

        s = '-'
        mut = (chro, pos, ref, alt)
        substitutionG = s.join(mut)

        try:
            qual = dict_qc[substitutionG]
        except Exception:
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

        # get list of EA for each gene
        if gene not in ea_dict:
            ea_dict[gene] = []
        if action in ['silent', '-', 'no_action', 'no_trace', 'no_gene'] or gene=='':
            continue
        elif action in ['STOP', 'no_STOP', 'STOP-loss', 'START_loss']:
            ea_dict[gene].append(1.)
        else:
            try:
                spl_EA = action.strip().split(';')
                new_spl_EA = delNoStr(spl_EA)
                if new_spl_EA == []:
                    continue
                else:
                    average = np.mean(new_spl_EA)
                    ea_dict[gene].append(average/100.)
            except AttributeError:
                ea_dict[gene].append(float(action)/100.)

    for gene_idx, gene_ in enumerate(total_gene_list):
        try:
            matrix_freq[pt_idx][gene_idx] = variant_counter
        except KeyError:
            matrix_freq[pt_idx][gene_idx] = 0

        try:
            ea_list = ea_dict[gene_]
            matrix_sum[pt_idx][gene_idx] = np.sum(ea_list)
        except KeyError:
            matrix_sum[pt_idx][gene_idx] = 0


if __name__ == '__main__':
    input_folder = '/Users/ywkim/rosinante/ADSP/iDEAL_input_folder/' # this would change depending on
                                                                     # where you have Rosinante mounted on
    output_folder = '/output/'  # for now the output files will be stored locally

    quality_file = input_folder + 'snvquality_detailed_jamie.csv'   # generated using the raw data from dbgap
    dict_qc = get_quality(quality_file)

    phenotype_file = input_folder + 'phs000572.v7.pht005179.v1.p4.c1.CaseControlEnrichedPhenotypesWES_y1.HMB-IRB.txt'
    ADe2, ADe3, ADe4, HCe2, HCe3, HCe4 = get_phenotype(phenotype_file)

    """
    Need to get the list of all the genes sequenced/observed
    """
    total_gene_set = set()
    trauma_folder = input_folder + 'after_QCfilter_jamie/'
    for filename in os.listdir(trauma_folder):
        pt_file = os.path.join(trauma_folder, filename)
        df = pd.read_csv(pt_file, sep='\t', index_col=None)
        genes = list(set(df['GENE'].values.tolist()))
        total_gene_set.update(genes)

    total_gene_list = list(sorted(total_gene_set))
    print(len(total_gene_list))

    group_names = ['ADe2', 'ADe3', 'ADe4', 'HCe2', 'HCe3', 'HCe4']
    for grp_idx, group in enumerate([ADe2, ADe3, ADe4, HCe2, HCe3, HCe4]): # each 'item' in this list is a list of patient (="group")
        matrix_freq = np.zeros((len(group), len(total_gene_list))) # patient by gene
        matrix_sum = np.zeros((len(group), len(total_gene_list)))

        for idx, pt_id in enumerate(group):
            filename = pt_id + '.trauma'
            pt_file = os.path.join(trauma_folder, filename)

            fill_matrix(pt_file, idx)

        matrix_out(matrix_freq, group, total_gene_list, output_folder, group_names[grp_idx] + 'count')
        matrix_out(matrix_sum, group, total_gene_list, output_folder, group_names[grp_idx] + '_sum')

