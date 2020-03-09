'''
Created on Jan 28, 2020

@author: ywkim
'''

import pandas as pd
from collections import Counter

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

        if race == '5' and eth == '0':

            if state == '0' or state == 0:
                if apoe == '44' or apoe == '34':
                    risk_list.append(sub_id)    # healthy with risk

            elif state == '1' or state == 1:
                if apoe == '22' or apoe == '23':    # AD with protection
                    prot_list.append(sub_id)
            else:
                pass


if __name__ == '__main__':
    input_folder = '/Users/ywkim/rosinante/ADSP/iDEAL_input_folder/' # this would change depending on
                                                                     # where you have Rosinante mounted on
    output_folder = ''  # for now the output files will be stored locally
    control_folder = input_folder + 'Control/'
    case_folder = input_folder + 'Case/'
    phenotype_file = input_folder + 'phs000572.v7.pht005179.v1.p4.c1.CaseControlEnrichedPhenotypesWES_y1.HMB-IRB.txt'

    quality_file = input_folder + 'snvquality_detailed_jamie.csv'   # generated using the raw data from dbgap
    dict_qc = get_quality(quality_file)

    print('-----Getting quality DONE-----')

    # Group patients into 'protective', 'risk'
    prot_list = []
    risk_list = []

    get_phenotype(phenotype_file)

    print(len(prot_list))
    print(len(risk_list))

    print('-----Getting phenotype DONE-----')

    total_gene_list = []

    totalVarRisk = Counter()
    totalVarWT = Counter()
    totalVarALL = Counter()

    EA_risk = {}
    EA_wt = {}
    EA_all = {}