"""
Created on June 4, 2020

@author: ywkim


"""

import pandas as pd
import numpy as np
import seaborn as sns
from scipy.stats import mstats
from matplotlib import pyplot as plt
from itertools import zip_longest
import os

from germline_analyses import get_list_from_csv



def delNoStr(currList):
    tempList = []
    for x in currList:
        try:
            tempList.append(float(x))
        except(ValueError):  # When There is a string
            continue
    return tempList


def processGermlineFile(gFile, EAList, candidate):
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

                else:
                    try:
                        spl_EA = action.strip().split(';')
                        new_spl_EA = delNoStr(spl_EA)
                        if new_spl_EA == []:
                            pass
                        else:
                            average = np.mean(new_spl_EA)
                            EAList.append(average)
                    except AttributeError:
                        EAList.append(float(action))


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

    if selection == 'two-sided':
        files.write('two-sided p-value is: ' + str(ks_P) + '\n')
    elif selection == 'less':
        files.write('positive selection p-value is: ' + str(ks_P) + '\n')
    else:
        files.write('negative selection p-value is: ' + str(ks_P) + '\n')


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


if __name__ == '__main__':
    phenotype_file = '/media/vision/ExtraDrive1/Exome/ADSP_discovery/' \
                     'phs000572.v7.pht005179.v1.p4.c1.CaseControlEnrichedPhenotypesWES_y1.HMB-IRB.txt'
    ADe2, ADe3, ADe4, HCe2, HCe3, HCe4 = get_phenotype(phenotype_file)

    cases = ADe2 + ADe3 + ADe4
    controls = HCe2 + HCe3 + HCe4

    trauma_folder = '/media/vision/ExtraDrive1/Exome/ADSP_discovery/all_short_name/'
    output_folder = '/home/vision/Documents/GermlineProject/ADSP/iDEAL_updated_test_2019/'

    EAListCase = []
    EAListControl = []

    # gene_file = output_folder + 'iDEAL_genelist.txt'
    # candidate = get_list_from_csv(gene_file, 'Gene', sep='\t')

    pathogenic_file = output_folder + 'pathogenic.txt'
    protective_file = output_folder + 'protective.txt'

    pathogenic_gene = get_list_from_csv(pathogenic_file, 'Gene', sep='\t')
    protective_gene = get_list_from_csv(protective_file, 'Gene', sep='\t')

    ptCounter = 0
    for filename in os.listdir(trauma_folder):
        ptFile = (os.path.join(trauma_folder, filename))
        if filename in cases:
            processGermlineFile(ptFile, EAListCase, pathogenic_gene)
        elif filename in controls:
            processGermlineFile(ptFile, EAListControl, pathogenic_gene)
        ptCounter += 1
        print(ptCounter)

    getHistogram(EAListCase, EAListControl, 'pathogenic_CaseControl')

    file1 = open(output_folder + 'pathogenic_case_control_kstest.txt', 'w')

    performKS_Test(EAListCase, EAListControl, 'less', file1)
    performKS_Test(EAListCase, EAListControl, 'greater', file1)
    performKS_Test(EAListCase, EAListControl, 'two-sided', file1)

    file1.close()


