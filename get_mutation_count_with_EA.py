"""
Created on 

@author: ywkim
"""

import os
import csv
import sys
import numpy as np
import pandas as pd
from collections import Counter

# from IPython import embed


def get_quality(qual_file):
    qual_dict = {}
    for line in open(qual_file):
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

        qual_dict[substitution] = filt
    return qual_dict


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
                    apoe3_hc.append(subID)
                if APOE == '34' or APOE =='44':
                    apoe4_hc.append(subID)  # healthy with risk

            elif state == '1' or state == 1:
                if APOE == '33':
                    apoe3_ad.append(subID)
                if APOE == '23' or APOE == '22':
                    apoe2_ad.append(subID)  # AD with protection
            else:
                pass


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

        try:
            qual = qc_dict[substitutionG]
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
    output_file = target_directory + '216_mutation_count_' + name + '.csv'
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

    ControlFolder = '/media/vision/ExtraDrive1/Exome/ADSP/Control'
    CaseFolder = '/media/vision/ExtraDrive1/Exome/ADSP/Case'
    phenoFile = '/media/vision/ExtraDrive1/Exome/ADSP/phs000572.v7.pht005179.v1.p4.c1.CaseControlEnrichedPhenotypesWES_y1.HMB-IRB.txt'

    qualFile = '/home/vision/Documents/GermlineProject/ADSP/snvquality_detailed_jamie.csv'
    qc_dict = get_quality(qualFile)

    print('-----Getting quality DONE-----')

    # Group patients into APOE2-AD, APOE4-HC, APOE3-AD, APOE4-AD
    apoe2_ad = []
    apoe4_hc = []
    apoe3_ad = []
    apoe3_hc = []

    getPheno(phenoFile)

    print('-----Getting phenotype DONE-----')

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

    # target_directory = '/home/vision/Documents/GermlineProject/ADSP/RVEA_BaylorPass_nonHisWhite_2v4_h5py_STARTLOSS100/'
    target_directory = '/home/vision/Documents/GermlineProject/ADSP/iDEAL_updated_test_2019/'
    gene_file = target_directory + 'iDEAL_genelist.txt'
    dfGene = pd.read_csv(gene_file, index_col=None)

    gene_list = dfGene['Gene'].values.tolist()

    print('-----Getting info-----')

    errorlog = []

    ptCount = 0
    for filename in os.listdir(ControlFolder):
        ptID = filename.split('.')[0]
        ptFile = (os.path.join(ControlFolder, filename))
        if ptID in apoe4_hc:
            getInfo(ptFile, ea_apoe4_hc, sub_apoe4_hc, het_count_apoe4_hc, homo_count_apoe4_hc)

        elif ptID in apoe3_hc:
            getInfo(ptFile, ea_apoe3_hc, sub_apoe3_hc, het_count_apoe3_hc, homo_count_apoe3_hc)

        else:
            pass
        ptCount += 1

    print(ptCount)

    for filename in os.listdir(CaseFolder):
        ptID = filename.split('.')[0]
        ptFile = (os.path.join(CaseFolder, filename))
        if ptID in apoe2_ad:
            getInfo(ptFile, ea_apoe2_ad, sub_apoe2_ad, het_count_apoe2_ad, homo_count_apoe2_ad)

        elif ptID in apoe3_ad:
            getInfo(ptFile, ea_apoe3_ad, sub_apoe3_ad, het_count_apoe3_ad, homo_count_apoe3_ad)

        else:
            pass
        ptCount += 1

    print(ptCount)

    output_dict(sub_apoe2_ad, het_count_apoe2_ad, homo_count_apoe2_ad, 'APOE2-AD')
    output_dict(sub_apoe4_hc, het_count_apoe4_hc, homo_count_apoe4_hc, 'APOE4-HC')
    output_dict(sub_apoe3_ad, het_count_apoe3_ad, homo_count_apoe3_ad, 'APOE3-AD')
    output_dict(sub_apoe3_hc, het_count_apoe3_hc, homo_count_apoe3_hc, 'APOE3-HC')

    # TODO This part needs to be re-written at some point
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

