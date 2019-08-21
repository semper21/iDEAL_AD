"""
Created on 

@author: ywkim
"""

#TODO Right now the comparison is between 3v3. Need to generalize this code to include 2v4 compariosn as well!

import os
import csv
import numpy as np
import pandas as pd
import seaborn as sns
from collections import Counter
import scipy.stats as stats
from matplotlib import pyplot as plt
import itertools as it

def getQuality(qualFile, dictQC):
    for line in open(qualFile):
        if line[0] == '#':
            continue
        cols = line.strip().split('\t')
        chr = cols[0]
        pos = cols[1]
        ref = cols[2]
        alt = cols[3]

        filt = cols[4]

        s = '-'
        mut = (chr, pos, ref, alt)
        substitution = s.join(mut)

        dictQC[substitution] = filt


def getPheno(phenoFile):
    for line in open(phenoFile):
        if line[0] == '#' or line[0] == 'd':
            continue
        cols = line.strip().split('\t')

        subID = cols[1]
        APOE = str(cols[7])
        race = str(cols[10])
        eth = str(cols[11])
        state = cols[13]
        # print type(APOE)

        if race == '5' and eth == '0':

            if state == '0' or state == 0:
                if APOE == '33':
                    Risklist.append(subID)  # healthy with risk


            elif state == '1' or state == 1:
                if APOE == '33':  # AD with protection
                    WTlist.append(subID)
            else:
                pass

    # print WTlist
    # print Risklist


def delNoStr(currList):
    tempList = []
    for x in currList:
        try:
            tempList.append(float(x))
        except(ValueError):  # When There is a string
            if (x in ['STOP', 'no_STOP', 'STOP-loss', 'START_loss']):
                tempList.append(1.)
            else:
                continue
    return tempList


def delNoStr_silent(currList):
    tempList = []
    for x in currList:
        try:
            tempList.append(float(x))
        except(ValueError):  # When There is a string
            if (x in ['STOP', 'no_STOP', 'STOP-loss', 'START_loss']):
                tempList.append(1.)
            elif (x in ['silent', '-', 'no_action', 'no_trace', 'no_gene']):
                tempList.append(0)
            else:
                print(x)
                continue
    return tempList


def getInfo(ptFile, EA_dict, count, EA_dict_silent, sub_counter):
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

        if gene == '':
            continue

        # get numbers for RVIS
        s = '-'
        mut = (chr, pos, ref, alt)
        substitutionG = s.join(mut)

        try:
            qual = dictQC[substitutionG]
        except Exception as e:
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

        # get sumEA
        if gene in geneList:
            # if gene == 'JMJD7':
            #    embed()
            count[gene] += 1
            if gene not in EA_dict:
                EA_dict[gene] = []
            if gene not in EA_dict_silent:
                EA_dict_silent[gene] = []

            if action in ['silent', '-', 'no_action', 'no_trace', 'no_gene']:
                EA_dict_silent[gene].append(0)

            elif action in ['STOP', 'no_STOP', 'STOP-loss', 'START_loss']:
                EA_dict[gene].append(1.)
                EA_dict_silent[gene].append(1.)
            else:
                try:
                    spl_EA = action.strip().split(';')
                    new_spl_EA = delNoStr(spl_EA)
                    new_spl_EA_silent = delNoStr_silent(spl_EA)
                    if new_spl_EA == []:
                        pass
                    else:
                        average = np.mean(new_spl_EA)
                        EA_dict[gene].append(average / 100.)

                    if new_spl_EA_silent == []:
                        pass
                    else:
                        average = np.mean(new_spl_EA_silent)
                        EA_dict_silent[gene].append(average / 100.)
                except AttributeError:
                    EA_dict[gene].append(float(action) / 100.)
                    EA_dict_silent[gene].append(float(action) / 100.)

            seq = sub + '(' + str(action) + ')'
            '''
            try:
                sub_counter[gene][sub] += 1
            except KeyError:
                sub_counter[gene] = Counter()
                sub_counter[gene][sub] += 1
            '''
            try:
                sub_counter[gene][seq] += 1
            except KeyError:
                sub_counter[gene] = Counter()
                sub_counter[gene][seq] += 1


def getBoxPlot(eacontrol, eacase, gene):
    EAListControl = list(np.asarray(eacontrol) * 100)
    EAListCase = list(np.asarray(eacase) * 100)

    avgEAcontrol = np.mean(EAListControl)
    avgEAcase = np.mean(EAListCase)

    EA = [EAListControl, EAListCase]

    rc = {'font.size': 20, 'axes.labelsize': 18, 'legend.fontsize': 16, 'axes.titlesize': 20, 'xtick.labelsize': 16,
          'ytick.labelsize': 16}
    sns.set(rc=rc)
    sns.set_style('whitegrid')
    ax = sns.violinplot(data=EA, palette='husl')
    # ax = sns.boxplot(data = EA, palette = 'husl')
    plt.ylim(0, 100)
    ax.set_xticklabels(['APOE2-AD', 'APOE4-HC'])
    plt.text(-0.1, 100, 'mean=' + '%.2f' % avgEAcontrol, fontsize=14)
    plt.text(0.9, 100, 'mean=' + '%.2f' % avgEAcase, fontsize=14)
    plt.title(gene)
    plt.savefig(targetDirc + gene + '_violin.png', dpi=200)
    plt.clf()


if __name__ == '__main__':

    ControlFolder = '/media/vision/ExtraDrive1/Exome/ADSP/Control'
    CaseFolder = '/media/vision/ExtraDrive1/Exome/ADSP/Case'
    phenoFile = '/media/vision/ExtraDrive1/Exome/ADSP/phs000572.v7.pht005179.v1.p4.c1.CaseControlEnrichedPhenotypesWES_y1.HMB-IRB.txt'

    qualFile = '/home/vision/Documents/GermlineProject/ADSP/snvquality_detailed_jamie.csv'
    dictQC = {}
    getQuality(qualFile, dictQC)

    print('-----Getting quality DONE-----')

    # Group patients into WT, Risk
    WTlist = []
    Risklist = []

    getPheno(phenoFile)

    print(len(WTlist))
    print(len(Risklist))

    print('-----Getting phenotype DONE-----')

    EA_Risk = {}
    EA_WT = {}
    count_Risk = Counter()
    count_WT = Counter()

    EA_Risk_silent = {}
    EA_WT_silent = {}

    targetDirc = '/home/vision/Documents/GermlineProject/ADSP/RVEA_BaylorPass_nonHisWhite_2v4_h5py_STARTLOSS100/'

    geneFile = targetDirc + 'iDEAL_genelist.txt'
    dfGene = pd.read_csv(geneFile, index_col=None)

    geneList = dfGene['Gene'].values.tolist()

    sub_risk = {}
    sub_WT = {}

    print('-----Getting info-----')

    errorlog = []

    ptCount = 0
    for filename in os.listdir(ControlFolder):
        ptID = filename.split('.')[0]
        ptFile = (os.path.join(ControlFolder, filename))
        if ptID in Risklist:
            getInfo(ptFile, EA_Risk, count_Risk, EA_Risk_silent, sub_risk)

        else:
            pass
        ptCount += 1
        print(ptCount)
    print(ptCount)

    for filename in os.listdir(CaseFolder):
        ptID = filename.split('.')[0]
        ptFile = (os.path.join(CaseFolder, filename))
        if ptID in WTlist:
            getInfo(ptFile, EA_WT, count_WT, EA_WT_silent, sub_WT)

        else:
            pass
        ptCount += 1
        print(ptCount)

    print(ptCount)
    '''
    outputFile = targetDirc + '216_genebygene_stats.csv'
    with open(outputFile, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(
            ['gene', 'case_mean', 'case_median', 'control_mean', 'control_median', 't-test', 'Welch', 'MWU test',
             'n_case', 'n_control'])

        for gene in geneList:

            try:
                EAlistRisk = EA_Risk[gene]
            except KeyError:
                EAlistRisk = [0]
            try:
                EAlistWT = EA_WT[gene]
            except KeyError:
                EAlistWT = [0]

            mean_case = np.mean(EAlistRisk)
            med_case = np.median(EAlistRisk)
            mean_control = np.mean(EAlistWT)
            med_control = np.median(EAlistWT)

            t = stats.ttest_ind(EAlistRisk, EAlistWT, equal_var=True)[1]
            w = stats.ttest_ind(EAlistRisk, EAlistWT, equal_var=False)[1]
            try:
                mwu = stats.mannwhitneyu(EAlistRisk, EAlistWT, use_continuity=False, alternative='two-sided')[1]
            except ValueError:
                mwu = '-'
            n = len(EAlistRisk)
            m = len(EAlistWT)
            info = ([gene, mean_case, med_case, mean_control, med_control, t, w, mwu, n, m])
            writer.writerow(info)
    '''

    dfRisk = pd.DataFrame.from_dict(sub_risk, orient='index')
    dfRisk = dfRisk.stack()
    dfRisk.to_csv(targetDirc + '216_mutation_count_APOE3-AD.txt', index=True, header=['Gene, Sub, Count'])

    dfWT = pd.DataFrame.from_dict(sub_WT, orient='index')
    dfWT = dfWT.stack()
    dfWT.to_csv(targetDirc + '216_mutation_count_APOE3-HC.txt', index=True, header=['Gene, Sub, Count'])

    # embed()


