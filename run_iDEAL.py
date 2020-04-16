"""
Created on Aug 14, 2019 (this new, 3.6 version)

@author: ywkim

-Started with old version (python 2.7)

"""

# TODO: need to clean this up (convert to python3)

import os
import csv
import sys
import numpy as np
import pandas as pd
import seaborn as sns
from random import shuffle
import matplotlib.pyplot as plt
from collections import Counter
from sklearn.linear_model import LinearRegression

import h5py


def getQuality(input_file):
    dict_ = {}
    for line in open(input_file):
        if line[0] == '#':
            continue
        cols = line.strip().split('\t')
        chro = cols[0]
        pos = cols[1]
        ref = cols[2]
        alt = cols[3]

        filt = cols[4]

        s = '-'
        mut = (chro, pos, ref, alt)
        substitution = s.join(mut)

        dict_[substitution] = filt

    return dict_


def getPhenotype(pheno_file):
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


def delNoStr(current_list):
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


def getInfo(ptFile, EA, totalVar, EA_ALL, totalVarALL, geneListALL):
    for line in open(ptFile):
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

        #get numbers for RVIS
        s = '-'
        mut = (chro, pos, ref, alt)
        substitutionG = s.join(mut)

        try:
            qual = dict_qc[substitutionG]
        except Exception as e:
            logf.write(filename + '\t' + str(e) + '\n')
            qual = 'non-PASS'

        if qual == 'non-PASS':
            continue
        elif qual == 'PASS':
            pass
        else:
            print('ERROR')

        if sub == '-': #no indels
            continue
        if gene == 'APOE' and sub == 'C130R':
            continue


        #get number of all variants (no maf cutoff, all var(except indels))
        if gene not in geneListALL:
            geneListALL.append(gene)
        totalVar[gene] += 1
        totalVarALL[gene] +=1

        #get sumEA
        if gene not in EA:
            EA[gene] = []
        if action in ['silent', '-', 'no_action', 'no_trace', 'no_gene'] or gene=='':
            continue
        elif action in ['STOP', 'no_STOP', 'STOP-loss', 'START_loss']:
            EA[gene].append(1.)
        else:
            try:
                spl_EA = action.strip().split(';')
                new_spl_EA = delNoStr(spl_EA)
                if new_spl_EA == []:
                    continue
                else:
                    average = np.mean(new_spl_EA)
                    EA[gene].append(average/100.)
            except AttributeError:
                EA[gene].append(float(action)/100.)

        if gene not in EA_ALL:
            EA_ALL[gene] = []
        if action in ['silent', '-', 'no_action', 'no_trace', 'no_gene'] or gene=='':
            continue
        elif action in ['STOP', 'no_STOP', 'STOP-loss', 'START_loss']:
            EA_ALL[gene].append(1.)
        else:
            try:
                spl_EA = action.strip().split(';')
                new_spl_EA = delNoStr(spl_EA)
                if new_spl_EA == []:
                    continue
                else:
                    average = np.mean(new_spl_EA)
                    EA_ALL[gene].append(average/100.)
            except AttributeError:
                EA_ALL[gene].append(float(action)/100.)


def getSumEA(anyDict):
    newDict={}
    remove=[]
    for gene in anyDict:
        eaList = anyDict[gene]
        if eaList == []:
            remove.append(gene)
        else:
            pass

        newDict[gene]=sum(eaList)

    for key in remove:
        del newDict[key]
    return newDict


def save2h5py(anyDict, name):
    f = h5py.File(target_dirc + name + '.hdf5', 'w')
    for gene, EA in anyDict.items():
        if type(EA) == int:
            f.create_dataset(gene, data=np.array([EA], dtype=np.float32))
        elif type(EA) == list:
            f.create_dataset(gene, data=np.array(EA, dtype=np.float32))
            #sys.exit()


def calcRVEA_OneFit(gene_list, totalVarALL, sumEA_ALL, totalVarRisk, sumEA_Risk, totalVarWT, sumEA_WT,  name):
    X_list = []
    Y_list = []
    for gene in gene_list:
        try:
            x = totalVarALL[gene]
        except KeyError:
            x = 0
        try:
            y = sumEA_ALL[gene]
        except KeyError:
            y = 0
        X_list.append(x)
        Y_list.append(y)

    # fit LR (y=ax+b)
    X_mat = np.transpose(np.matrix(X_list))
    Y_mat = np.transpose(np.matrix(Y_list))
    linreg = LinearRegression()     # no need to fit to origin

    linreg.fit(X_mat, Y_mat)
    a = linreg.coef_[0][0]
    c = linreg.intercept_[0]

    print('a=',a, 'c=',c)

    RVEA_Risk = {}
    X_Risk = []
    Y_Risk = []
    outputFileAA = target_dirc + 'RVEA_Risk' + name
    with open(outputFileAA, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Gene', '#allMut', 'sumEA', 'Residual'])
        # ax+by+c=0 (b=-1)
        # r = ax-y+c/sqrt(a^2+b^2)
        for gene in gene_list:
            try:
                x = totalVarRisk[gene]
            except KeyError:
                x = 0

            try: #######################
                y = sumEA_Risk[gene]
            except KeyError:
                y = 0

            y_hat = (a*x + c)
            r = y - y_hat
            RVEA_Risk[gene] = r

            X_Risk.append(x)
            Y_Risk.append(y)


            info = ([gene, x, y, r])
            writer.writerow(info)

    getRegPlot(X_list, Y_list, X_Risk, Y_Risk, '# of all observed variants', 'EAburden', 'RVEA_Risk', name)


    RVEA_WT = {}
    X_WT = []
    Y_WT = []
    outputFileGG = target_dirc + 'RVEA_WT' + name
    with open(outputFileGG, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Gene', '#allMut', 'sumEA', 'Residual'])
        # ax+by+c=0 (b=-1)
        # r = ax-y+c/sqrt(a^2+b^2)
        for gene in gene_list:
            try:
                x = totalVarWT[gene]
            except KeyError:
                x = 0

            try: #######################
                y = sumEA_WT[gene]
            except KeyError:
                y = 0

            y_hat = (a*x + c)
            r = y - y_hat
            RVEA_WT[gene] = r

            X_WT.append(x)
            Y_WT.append(y)

            info = ([gene, x, y, r])
            writer.writerow(info)

    getRegPlot(X_list, Y_list, X_WT, Y_WT, '# of all observed variants', 'EAburden', 'RVEA_WT', name)

    return RVEA_Risk, RVEA_WT


def calcRofRVEA(gene_list, r_risk, r_wt, name): # this is for only genes with mutations in both cohorts
    outputFile = target_dirc + 'ADSP_RofRVEA' + name
    ideal = []
    dict_residual = {}
    with open(outputFile, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Gene', 'risk_residual', 'wt_residual', 'iDEAL'])
        X_list = []
        Y_list = []
        for gene in gene_list:

            x = r_risk[gene]
            y = r_wt[gene]

            X_list.append(x)
            Y_list.append(y)

        X_mat = np.transpose(np.matrix(X_list))
        Y_mat = np.transpose(np.matrix(Y_list))
        linreg = LinearRegression()

        linreg.fit(X_mat, Y_mat)
        a = linreg.coef_[0][0]
        c = linreg.intercept_[0]

        print('a = ', a, 'c = ', c)
        # ax+by+c=0 (b=-1)
        # r = ax-y+c/sqrt(a^2+b^2)
        for gene in gene_list:
            x = r_risk[gene]
            y = r_wt[gene]

            y_hat = (a*x + c)
            r = y - y_hat

            if name == '_real':
                dict_residual[gene] = r
            elif name == '_random':
                if gene not in dict_residual:
                    dict_residual[gene] = []
                dict_residual[gene].append(r)
            else:
                print(name)
                ideal.append(r)
            info = ([gene, x, y, r])
            writer.writerow(info)

    getRegPlot(X_list, Y_list, X_list, Y_list, 'RVEA_Risk', 'RVEA_WT', 'residualOfRVEAs', name)
    return dict_residual, ideal


def getRegPlot(x_reg, y_reg, xx, yy, x_label, y_label, title, name):
    plt.axhline(0, color='0.75')
    plt.axvline(0, color='0.75')
    sns.regplot(x = np.asarray(xx), y=np.asarray(yy), fit_reg=False)
    ax = sns.regplot(x = np.asarray(x_reg), y=np.asarray(y_reg), scatter=False)
    ax.set(xlabel=x_label, ylabel=y_label)
    if y_label == 'EAburden':
        plt.xlim(xmin=0)
        plt.ylim(ymin=0)
    else:
        pass
    plt.title(title)
    plt.savefig(target_dirc + title + name + '.png', dpi = 200)
    plt.clf()


if __name__ == '__main__':
    control_folder = '/media/vision/ExtraDrive1/Exome/ALZ/ADSP/Control'
    case_folder = '/media/vision/ExtraDrive1/Exome/ALZ/ADSP/Case'
    phenotype_file = '/media/vision/ExtraDrive1/Exome/ALZ/ADSP/phs000572.v7.pht005179.v1.p4.c1.' \
                     'CaseControlEnrichedPhenotypesWES_y1.HMB-IRB.txt'

    quality_file = '/home/vision/Documents/GermlineProject/ADSP/snvquality_detailed_jamie.csv'
    dict_qc = getQuality(quality_file)

    print('-----Getting quality DONE-----')

    # Group patients into WT, Risk
    prot_list = []
    risk_list = []

    getPhenotype(phenotype_file)

    print(len(prot_list))
    print(len(risk_list))

    print('-----Getting phenotype DONE-----')

    total_gene_list = []

    # h5py
    totalVarRisk = Counter()
    totalVarWT = Counter()
    totalVarALL = Counter()

    EA_risk = {}
    EA_wt = {}
    EA_all = {}

    target_dirc = '/home/vision/Documents/GermlineProject/ADSP/RVEA_BaylorPass_nonHisWhite_2v4_h5py_STARTLOSS100/'

    print('-----Getting info-----')

    ptCount = 0
    outputFile = target_dirc + 'RiskAlleleStatus.csv'
    logf = open(target_dirc + 'error.log', 'w')
    with open(outputFile, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['PT_ID', 'AlleleStatus'])
        for filename in os.listdir(control_folder):
            ptID = filename.split('.')[0]
            ptFile = (os.path.join(control_folder, filename))
            if ptID in risk_list:
                status = 'Risk'
                getInfo(ptFile, EA_risk, totalVarRisk, EA_all, totalVarALL, total_gene_list)

            else:
                status = 'Healthy_Something_else'
            ptCount += 1
            print(ptCount)
            info = ([ptID, status])
            writer.writerow(info)
        print(ptCount)

        for filename in os.listdir(case_folder):
            ptID = filename.split('.')[0]
            ptFile = (os.path.join(case_folder, filename))
            if ptID in prot_list:
                status = 'WT'
                getInfo(ptFile, EA_wt, totalVarWT, EA_all, totalVarALL, total_gene_list)

            else:
                status = 'Case_Something_else'

            ptCount += 1
            print(ptCount)
            info = ([ptID, status])
            writer.writerow(info)

        print(ptCount)

    # h5py
    sumEA_Risk = getSumEA(EA_risk)
    sumEA_WT = getSumEA(EA_wt)
    sumEA_ALL = getSumEA(EA_all)

    save2h5py(totalVarRisk, 'totalVarRisk')
    save2h5py(totalVarWT, 'totalVarWT')
    save2h5py(totalVarALL, 'totalVarALL')
    save2h5py(sumEA_Risk, 'sumEA_Risk')
    save2h5py(sumEA_WT, 'sumEA_WT')
    save2h5py(sumEA_ALL, 'sumEA_ALL')

    print('-----Initializing RVEA calculation-----')
    # calculate RVEA of AA and GG using same fitted line (of all GG, AA, GA people)
    normRVEA_Risk, normRVEA_WT = calcRVEA_OneFit(total_gene_list, totalVarALL, sumEA_ALL, totalVarRisk, sumEA_Risk,
                                                 totalVarWT, sumEA_WT, '_real')
    save2h5py(normRVEA_Risk, 'RVEA_Risk')
    save2h5py(normRVEA_WT, 'RVEA_WT')

    RofRVEA, RVEA_list = calcRofRVEA(total_gene_list, normRVEA_Risk, normRVEA_WT, '_real')

    save2h5py(RofRVEA, 'RofRVEA')

    print('-----Initializing random RVEA calculation-----')
    ################
    #    RANDOM    #
    ################

    totalPt = risk_list + prot_list

    randAAlist_total = []
    randGGlist_total = []
    for p in range(100):
        randAAlist = []
        randGGlist = []
        shuffle(totalPt)

        for i in range(len(totalPt)):
            pt = totalPt[i]  # pt = each patient
            # ptN = germlineFiles[r].rsplit('/',1)[1]
            # pt = ptN.split('.')[0]
            if i < len(risk_list):  # Risk
                randAAlist.append(pt)
            else:
                randGGlist.append(pt)

        randAAlist_total.append(randAAlist)
        randGGlist_total.append(randGGlist)

    df = pd.DataFrame(randAAlist_total)
    df.to_csv(target_dirc + 'randomptsetsRisk.csv', sep='\t', index=False)

    df = pd.DataFrame(randGGlist_total)
    df.to_csv(target_dirc + 'randomptsetsWT.csv', sep='\t', index=False)

    RVISofRVIS_random = {}
    for ptidx in range(100):

        RgeneListTotal = []

        RtotalVarAA = Counter()
        RtotalVarGG = Counter()
        RtotalVarALL = Counter()

        RandEA_AA = {}
        RandEA_GG = {}
        RandEA_ALL = {}

        perc = 0

        randAAlist = randAAlist_total[ptidx]
        randGGlist = randGGlist_total[ptidx]

        for pt in totalPt:
            if pt in risk_list:
                ptFile = (os.path.join(control_folder, pt))
            elif pt in prot_list:
                ptFile = (os.path.join(case_folder, pt))
            else:
                print('error')
                sys.exit()

            if pt in randAAlist:
                status = 'Risk'
                getInfo(ptFile, RandEA_AA, RtotalVarAA, RandEA_ALL, RtotalVarALL, RgeneListTotal)

            elif pt in randGGlist:
                status = 'WT'
                getInfo(ptFile, RandEA_GG, RtotalVarGG, RandEA_ALL, RtotalVarALL, RgeneListTotal)

            else:
                status = 'GA'
            perc += 1
            if perc == 1000:
                print('1000 random people done')

        RsumEA_AA = getSumEA(RandEA_AA)
        RsumEA_GG = getSumEA(RandEA_GG)
        RsumEA_ALL = getSumEA(RandEA_ALL)

        noMutGene = []
        # calculate RVIS of AA and GG using same fitted line (of all GG, AA, GA people)
        normRVIS_AA_random, normRVIS_GG_random = calcRVEA_OneFit(RgeneListTotal, RtotalVarALL, RsumEA_ALL, RtotalVarAA,
                                                                 RsumEA_AA, RtotalVarGG, RsumEA_GG, '_random')
        # uses 'normalized' RVIS (one fitted line for calculating RVIS of both AA and GG)
        RVISofRVIS_random, RVIS_list_random = calcRofRVEA(RgeneListTotal, normRVIS_AA_random, normRVIS_GG_random,
                                                          RVISofRVIS_random, '_random')

        print('random # = ', ptidx)

    save2h5py(RVISofRVIS_random, 'RofRVEA_random')

    # Get z-scores
    noInfo = 0
    outputFile = target_dirc + 'ControlledZscores'
    with open(outputFile, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Gene', 'residual', 'mean', 'std', 'z-score', 'RofRVEA_random'])
        for gene in total_gene_list:
            try:
                mean = np.mean(RVISofRVIS_random[gene])
                std = np.std(RVISofRVIS_random[gene])
                x = RofRVEA[gene]
                z = (x - mean) / std
                r = RVISofRVIS_random[gene]
                new_r = []
                # for d in r:
                #    new_r.append(float('{0:.2f}'.format(d)))
                # info = ([gene, float('{0:.2f}'.format(x)), float('{0:.2f}'.format(mean)),
                #          float('{0:.2f}'.format(std)), float('{0:.2f}'.format(z)), new_r])

                info = ([gene, x, mean, std, z, r])
                writer.writerow(info)
            except KeyError:
                noInfo += 1

    print('number of genes with no z-score = ', noInfo)
