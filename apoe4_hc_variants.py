"""
Created on Aug 13, 2020

@author: ywkim
"""

import os
import csv
from pathlib import Path
from collections import Counter
from germline_analyses import get_list_from_csv


def get_phenotype(pheno_file):
    AD2_hom, AD2_het, AD3, AD4_hom, AD4_het = [], [], [], [], []
    HC2_hom, HC2_het, HC3, HC4_hom, HC4_het = [], [], [], [], []

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
                if apoe == '22':
                    AD2_hom.append(sub_id)
                elif apoe == '23':
                    AD2_het.append(sub_id)
                elif apoe == '33':
                    AD3.append(sub_id)
                elif apoe == '44':
                    AD4_hom.append(sub_id)  # healthy with risk
                elif apoe == '34':
                    AD4_het.append(sub_id)
                else:
                    pass

            elif str(state) == '0':  # HC
                if apoe == '22':
                    HC2_hom.append(sub_id)
                elif apoe == '23':
                    HC2_het.append(sub_id)
                elif apoe == '33':  # APOE3
                    HC3.append(sub_id)
                elif apoe == '44':  # APOE4
                    HC4_hom.append(sub_id)
                elif  apoe == '34':
                    HC4_het.append(sub_id)

                else:
                    pass
            else:
                pass

    # return AD2_hom, AD2_het, AD3, AD4_hom, AD4_het, HC2_hom, HC2_het, HC3, HC4_hom, HC4_het
    return AD2_hom, HC2_hom, HC2_het, HC4_hom, AD4_hom, AD4_het


def delNoStr(currList):
    tempList = []
    for x in currList:
        try:
            tempList.append(float(x))
        except(ValueError):  # When There is a string
            continue
    return tempList


def extract_variants(file_, dict_, genelist_):
    for line in open(file_):
        if line[0] == '#':
            continue
        cols = line.strip().split('\t')
        gene = cols[4].split(';')[0]
        sub = cols[5]
        action = cols[6]

        if gene == '':
            continue

        if sub == '-':  # no indels
            continue

        if gene in genelist_:
            if action in ['silent', '-', 'no_action', 'no_trace', 'no_gene', 'indel', 'fs-indel'] or gene == '':
                continue

            else:
                if action in ['STOP', 'no_STOP', 'STOP-loss', 'START_loss']:
                    try:
                        dict_[gene].append(sub)
                    except KeyError:
                        dict_[gene] = []
                        dict_[gene].append(sub)
                else:
                    try:
                        spl_EA = action.strip().split(';')
                        new_spl_EA = delNoStr(spl_EA)
                        if new_spl_EA == []:
                            pass
                        else:
                            try:
                                dict_[gene].append(sub)
                            except KeyError:
                                dict_[gene] = []
                                dict_[gene].append(sub)
                    except AttributeError:
                        try:
                            dict_[gene].append(sub)
                        except KeyError:
                            dict_[gene] = []
                            dict_[gene].append(sub)


if __name__ == '__main__':
    input_folder = str(Path().absolute()) + '/input/'
    output_folder = str(Path().absolute()) + '/output/'

    # gene_file = input_folder + 'iDEAL_genelist.txt'
    # gene_file = input_folder + 'top_pathogenic.txt'
    gene_file = input_folder + 'top_protective.txt'
    gene_list = get_list_from_csv(gene_file, 'Gene', sep='\t')

    phenotype_file = '/media/vision/ExtraDrive1/Exome/ADSP_discovery/' \
                     'phs000572.v7.pht005179.v1.p4.c1.CaseControlEnrichedPhenotypesWES_y1.HMB-IRB.txt'

    patient_lists = get_phenotype(phenotype_file)[-3:]
    # labels = ['AD2_hom', 'HC2_hom', 'HC2_het']
    labels = ['HC4_hom', 'AD4_hom', 'AD4_het']

    trauma_folder = '/media/vision/ExtraDrive1/Exome/ADSP_discovery/all_short_name/'

    for idx, patient_list in enumerate(patient_lists):
        variant_dict = {}
        ptCounter = 0
        label = labels[idx]
        for filename in os.listdir(trauma_folder):
            ptFile = (os.path.join(trauma_folder, filename))
            if filename in patient_list:
                extract_variants(ptFile, variant_dict, gene_list)
            ptCounter += 1
            print(ptCounter)


        output_file = output_folder + label + '_substitution_list.txt'

        with open(output_file, 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['Gene', 'Sub', 'Count'])
            flag = []
            for gene in gene_list:
                try:
                    sub_counter = Counter(variant_dict[gene])
                    for sub in sub_counter:
                        info = [gene, sub, sub_counter[sub]]
                        writer.writerow(info)
                except KeyError:
                    pass

