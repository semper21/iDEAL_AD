"""
Created on Aug 3, 2020

@author: ywkim

Trying to compare mutational burden in "protective" iDEAL genes in healthy e4 homozygotes vs heterozygotes
"""

import os
import csv
import pandas as pd
import numpy as np
import seaborn as sns
from scipy.stats import mstats
from matplotlib import pyplot as plt

from pathlib import Path
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
    return HC4_hom, HC4_het


def delNoStr(currList):
    tempList = []
    for x in currList:
        try:
            tempList.append(float(x))
        except(ValueError):  # When There is a string
            continue
    return tempList


def processGermlineFile(gFile, EAList, ea_dict, candidate):
    for line in open(gFile):
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

        if gene in candidate:
            if action in ['silent', '-', 'no_action', 'no_trace', 'no_gene', 'indel', 'fs-indel'] or gene == '':
                continue

            else:
                if action in ['STOP', 'no_STOP', 'STOP-loss', 'START_loss']:
                    EAList.append(100.0)
                    try:
                        ea_dict[gene].append(100.0)
                    except KeyError:
                        ea_dict[gene] = []
                        ea_dict[gene].append(100.0)
                else:
                    try:
                        spl_EA = action.strip().split(';')
                        new_spl_EA = delNoStr(spl_EA)
                        if new_spl_EA == []:
                            pass
                        else:
                            average = np.mean(new_spl_EA)
                            EAList.append(average)
                            try:
                                ea_dict[gene].append(average)
                            except KeyError:
                                ea_dict[gene] = []
                                ea_dict[gene].append(average)
                    except AttributeError:
                        EAList.append(float(action))
                        try:
                            ea_dict[gene].append(float(action))
                        except KeyError:
                            ea_dict[gene] = []
                            ea_dict[gene].append(float(action))



def getFraction(EAD, length):
    fracList = []
    for ea in EAD:
        f = float(ea) / length
        fracList.append(f)
    return fracList


def plotHistogram(EADCase, EADControl, name, types):
    bin = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

    df = pd.DataFrame({'Case': EADCase, 'Control': EADControl})
    df['bin'] = bin
    df = pd.melt(df, id_vars="bin", value_vars=['Case', 'Control'], var_name='State', value_name='Count')
    # print df

    rc = {'font.size': 12, 'axes.labelsize': 16, 'legend.fontsize': 10, 'axes.titlesize': 20, 'xtick.labelsize': 14,
          'ytick.labelsize': 14}
    sns.set(rc=rc)
    sns.set_style('whitegrid')
    # plt.figure(figsize = (15, 10))

    ax = sns.barplot(x='bin', y='Count', hue='State', data=df, palette='husl')
    ax.set(xlabel='EA Scores', ylabel=types)
    ax.set_xticklabels([10, 20, 30, 40, 50, 60, 70, 80, 90, 100], rotation=0)

    plt.title('EAD of Germline Variants (' + name + ')')
    plt.savefig(output_folder + name + '_' + types + '.png', dpi=300)
    plt.clf()


def getHistogram(caseList, controlList, name):
    EADCase = np.histogram(caseList, bins=[0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
    EADCase = [item.tolist() for item in EADCase]
    EADCase = EADCase[0]

    EADControl = np.histogram(controlList, bins=[0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
    EADControl = [item.tolist() for item in EADControl]
    EADControl = EADControl[0]

    EADCaseF = getFraction(EADCase, len(caseList))
    EADControlF = getFraction(EADControl, len(controlList))

    plotHistogram(EADCase, EADControl, name, 'Frequency')
    plotHistogram(EADCaseF, EADControlF, name, 'Fraction')

    bin = [10,20,30,40,50,60,70,80,90,100]

    df = pd.DataFrame({'Bin' : bin, 'Case': EADCase, 'Control': EADControl})
    df.to_csv(output_folder + 'EAD_bins' + name + '.csv')

def performKS_Test(caseList, controlList, selection, files):
    ks_P = mstats.ks_twosamp(caseList, controlList, alternative=selection)[1]

    if files:
        if selection == 'two-sided':
            files.write('two-sided p-value is: ' + str(ks_P) + '\n')
        elif selection == 'less':
            files.write('positive selection p-value is: ' + str(ks_P) + '\n')
        else:
            files.write('negative selection p-value is: ' + str(ks_P) + '\n')

    return ks_P


if __name__ == '__main__':
    input_folder = str(Path().absolute()) + '/input/'
    output_folder = str(Path().absolute()) + '/output_ADSP_discovery/'

    # gene_file = input_folder + 'protective.txt'
    # gene_file = input_folder + 'iDEAL_genelist.txt'
    # gene_list = get_list_from_csv(gene_file, 'Gene', sep='\t')

    gene_list = ['BMS1', 'PLXNA4', 'STEAP1B', 'VARS', 'ATP1A1', 'MEN1', 'B3GAT2', 'ARHGAP33', 'IGFALS', 'RTN1', 'ZNF20',
                 'TTC7B', 'MAN1A1', 'KIF5A', 'MAST1', 'WSB2', 'SSTR3', 'LAPTM4B', 'NRXN2', 'SOX6', 'PTBP1', 'KIF3C',
                 'SATB1', 'GCLC', 'CYP27B1', 'ATXN1', 'NDRG2', 'PAK2', 'LAMA2', 'PANX2', 'SLC35B4', 'LRRC17', 'MASP1',
                 'BCL2L13', 'POMT1', 'CENPT', 'SMTNL1', 'SCAF11', 'RB1CC1', 'PRPH', 'PTCH1', 'TRAF3IP2', 'OR52E6',
                 'SYTL2', 'STARD7', 'UGT3A1', 'DDIT4L', 'SLCO3A1', 'ALKBH4', 'EPHA7', 'GPIHBP1', 'CRTC3', 'OPRD1',
                 'YDJC', 'COX11', 'TBX1', 'ZCCHC9', 'PELO', 'ZAR1', 'TUBB3', 'ABHD2', 'ANAPC15', 'TMEM163', 'RBM17',
                 'SLC35G6', 'CMC1', 'SMARCD2', 'HIST1H3G']

    phenotype_file = '/media/vision/ExtraDrive1/Exome/ADSP_discovery/' \
                     'phs000572.v7.pht005179.v1.p4.c1.CaseControlEnrichedPhenotypesWES_y1.HMB-IRB.txt'

    trauma_folder = '/media/vision/ExtraDrive1/Exome/ADSP_discovery/all_short_name/'

    EAListCase = []
    EAListControl = []
    EA_dict_case = {}
    EA_dict_control = {}

    HC4_hom, HC4_het = get_phenotype(phenotype_file)

    cases = HC4_hom
    controls = HC4_het

    ptCounter = 0
    for filename in os.listdir(trauma_folder):
        ptFile = (os.path.join(trauma_folder, filename))
        if filename in cases:
            processGermlineFile(ptFile, EAListCase, EA_dict_case, gene_list)
        elif filename in controls:
            processGermlineFile(ptFile, EAListControl, EA_dict_control, gene_list)
        ptCounter += 1
        print(ptCounter)

    """
    # all genes
    getHistogram(EAListCase, EAListControl, 'APOE4_HC_hom_vs_het_all')

    file1 = open(output_folder + 'APOE4_HC_hom_vs_het_protective_kstest.txt', 'w')

    performKS_Test(EAListCase, EAListControl, 'less', file1)
    performKS_Test(EAListCase, EAListControl, 'greater', file1)
    performKS_Test(EAListCase, EAListControl, 'two-sided', file1)

    """
    # gene-level
    output_file = output_folder + 'APOE4_HC_hom_vs_het_protective_gene_by_gene_kstest.txt'

    with open(output_file, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Gene', 'test', 'p-value', 'hom', 'het'])
        for gene in gene_list:
            try:
                case_list = EA_dict_case[gene]
            except KeyError:
                print(gene, 'no_hom')
                case_list = [0]
            try:
                control_list = EA_dict_control[gene]
            except KeyError:
                control_list = 0
                print(gene, 'no_het')

            p_less = performKS_Test(case_list, control_list, 'less', None)
            # p_greater = performKS_Test(case_list, control_list, 'greater', None)
            # p_two = performKS_Test(case_list, control_list, 'two-sided', None)

            writer.writerow([gene, 'less', p_less, case_list, control_list])
            # writer.writerow([gene, 'greater', p_greater])
            # writer.writerow([gene, 'two_sided', p_two])

    """
    output_file = output_folder + 'APOE4_HC_hom_vs_het_protective_iterative_kstest.txt'

    gene_list_r = gene_list
    with open(output_file, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Gene_list', 'test', 'p-value', 'hc_e4_hom', 'hc_e4_het'])
        for i in range(len(gene_list_r)):
            case_list = []
            control_list = []
            for gene in gene_list_r:
                try:
                    case_list += EA_dict_case[gene]
                except KeyError:
                    case_list += [0]
                try:
                    control_list += EA_dict_control[gene]
                except KeyError:
                    control_list += [0]

            p_greater = performKS_Test(case_list, control_list, 'greater', None)
            # p_two = performKS_Test(case_list, control_list, 'two-sided', None)
            writer.writerow([gene_list_r, 'greater', p_greater, case_list, control_list])
            # writer.writerow([gene_list_r, 'two_sided', p_two])

            del gene_list_r[-1]


    """




