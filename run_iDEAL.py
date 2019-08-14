"""
Created on 

@author: ywkim

-Started with old version (python 2.7)

"""


import os
import csv
import glob
import sys
import numpy as np
import pandas as pd
import seaborn as sns
from random import shuffle
import matplotlib.pyplot as plt
from collections import Counter
from sklearn.linear_model import LinearRegression

import h5py

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
    for line in file(phenoFile):
        if line[0] == '#' or line[0] == 'd':
            continue
        cols = line.strip().split('\t')

        subID = cols[1]
        APOE = str(cols[7])
        race = str(cols[10])
        eth = str(cols[11])
        state = cols[13]
        #print type(APOE)

        if race == '5' and eth == '0':

            if state == '0' or state == 0:
                if APOE == '44' or APOE == '34':
                    Risklist.append(subID) #healthy with risk


            elif state == '1' or state == 1:
                if APOE == '22' or APOE == '23': #AD with protection
                    WTlist.append(subID)
            else:
                pass

    #print WTlist
    #print Risklist


def delNoStr(currList):
    tempList = []
    for x in currList:
        try:
            tempList.append(float(x))
        except(ValueError):  # When There is a string
            if(x in ['STOP', 'no_STOP', 'STOP-loss', 'START_loss']):
                tempList.append(1.0)
            else:
                continue
    return tempList

def getInfo(ptFile, EA, totalVar, EA_ALL, totalVarALL, geneListALL):
    for line in file(ptFile):
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

        #get numbers for RVIS
        s = '-'
        mut = (chr, pos, ref, alt)
        substitutionG = s.join(mut)

        try:
            qual = dictQC[substitutionG]
        except Exception as e:
            #error = (filename + '\t' + str(e))
            #errorlog.append(error)
            logf.write(filename + '\t' + str(e) + '\n')
            qual = 'non-PASS'

        if qual == 'non-PASS':
            continue
        elif qual == 'PASS':
            pass
        else:
            print 'ERROR'

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
    f = h5py.File(targetDirc + name + '.hdf5', 'w')
    for gene, EA in anyDict.items():
        if type(EA) == int:
            f.create_dataset(gene, data=np.array([EA], dtype=np.float32))
        elif type(EA) == list:
            f.create_dataset(gene, data=np.array(EA, dtype=np.float32))
            #sys.exit()

def calcRVEA_OneFit(geneList, totalVarALL, sumEA_ALL, totalVarRisk, sumEA_Risk, totalVarWT, sumEA_WT,  name):

    X_list = []
    Y_list = []
    for gene in geneList:
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

    #fit LR (y=ax+b)
    X_mat = np.transpose(np.matrix(X_list))
    Y_mat = np.transpose(np.matrix(Y_list))
    linreg = LinearRegression() #no need to fit to origin

    linreg.fit(X_mat, Y_mat)
    a = linreg.coef_[0][0]
    c = linreg.intercept_[0]

    print 'a=',a, 'c=',c

    RVEA_Risk = {}
    X_Risk = []
    Y_Risk = []
    outputFileAA = targetDirc + 'RVEA_Risk' + name
    with open(outputFileAA, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Gene', '#allMut', 'sumEA', 'Residual'])
        #ax+by+c=0 (b=-1)
        #r = ax-y+c/sqrt(a^2+b^2)
        for gene in geneList:
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
    outputFileGG = targetDirc + 'RVEA_WT' + name
    with open(outputFileGG, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Gene', '#allMut', 'sumEA', 'Residual'])
        #ax+by+c=0 (b=-1)
        #r = ax-y+c/sqrt(a^2+b^2)
        for gene in geneList:
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

def calcRofRVEA(geneList, RVEA_Risk, RVEA_WT, RofRVEA, name): #this is for only genes with mutations in both AA and GG
    outputFile = targetDirc + 'ADSP_RofRVEA' + name
    RVIS = []
    with open(outputFile, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Gene', 'RVIS of Risk', 'RVIS of WT', 'RVIS of RVIS'])
        X_list = []
        Y_list = []
        for gene in geneList:

            x = RVEA_Risk[gene]
            y = RVEA_WT[gene]

            X_list.append(x)
            Y_list.append(y)

        X_mat = np.transpose(np.matrix(X_list))
        Y_mat = np.transpose(np.matrix(Y_list))
        linreg = LinearRegression()

        linreg.fit(X_mat, Y_mat)
        a = linreg.coef_[0][0]
        c = linreg.intercept_[0]

        print 'a = ', a, 'c = ', c
        #ax+by+c=0 (b=-1)
        #r = ax-y+c/sqrt(a^2+b^2)
        for gene in geneList:
            x = RVEA_Risk[gene]
            y = RVEA_WT[gene]

            y_hat = (a*x + c)
            r = y - y_hat

            if name == '_real':
                RofRVEA[gene] = r
            elif name == '_random':
                if gene not in RofRVEA:
                    RofRVEA[gene] = []
                RofRVEA[gene].append(r)
            else:
                print name
            RVIS.append(r)
            info = ([gene, x, y, r])
            writer.writerow(info)

    getRegPlot(X_list, Y_list, X_list, Y_list, 'RVEA_Risk', 'RVEA_WT', 'residualOfRVEAs', name)
    return RofRVEA, RVIS

def getRegPlot(x_reg, y_reg, xx, yy, xlabel, ylabel, title, name):
    plt.axhline(0, color = '0.75')
    plt.axvline(0, color = '0.75')
    sns.regplot(x = np.asarray(xx), y = np.asarray(yy), fit_reg = False)
    ax = sns.regplot(x = np.asarray(x_reg), y = np.asarray(y_reg), scatter = False)
    ax.set(xlabel=xlabel, ylabel=ylabel)
    if ylabel == 'EAburden':
        plt.xlim(xmin=0)
        plt.ylim(ymin=0)
    else:
        pass
    plt.title(title)
    plt.savefig(targetDirc + title + name + '.png', dpi = 200)
    plt.clf()

if __name__ == '__main__':
    ControlFolder = '/media/vision/ExtraDrive1/Exome/ALZ/ADSP/Control'
    CaseFolder = '/media/vision/ExtraDrive1/Exome/ALZ/ADSP/Case'
    phenoFile = '/media/vision/ExtraDrive1/Exome/ALZ/ADSP/phs000572.v7.pht005179.v1.p4.c1.CaseControlEnrichedPhenotypesWES_y1.HMB-IRB.txt'

    qualFile = '/home/vision/Documents/GermlineProject/ADSP/snvquality_detailed_jamie.csv'
    dictQC = {}
    getQuality(qualFile, dictQC)

    print
    '-----Getting quality DONE-----'

    # Group patients into WT, Risk
    WTlist = []
    Risklist = []

    getPheno(phenoFile)

    print
    len(WTlist)
    print
    len(Risklist)

    print
    '-----Getting phenotype DONE-----'

    geneListTotal = []

    ##h5py
    totalVarRisk = Counter()
    totalVarWT = Counter()
    totalVarALL = Counter()

    EA_Risk = {}
    EA_WT = {}
    EA_ALL = {}

    targetDirc = '/home/vision/Documents/GermlineProject/ADSP/RVEA_BaylorPass_nonHisWhite_2v4_h5py_STARTLOSS100/'

    print
    '-----Getting info-----'

    errorlog = []

    ptCount = 0
    outputFile = targetDirc + 'RiskAlleleStatus.csv'
    logf = open(targetDirc + 'error.log', 'w')
    with open(outputFile, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['PT_ID', 'AlleleStatus'])
        for filename in os.listdir(ControlFolder):
            ptID = filename.split('.')[0]
            ptFile = (os.path.join(ControlFolder, filename))
            if ptID in Risklist:
                status = 'Risk'
                getInfo(ptFile, EA_Risk, totalVarRisk, EA_ALL, totalVarALL, geneListTotal)

            else:
                status = 'Healthy_Something_else'
            ptCount += 1
            print
            ptCount
            info = ([ptID, status])
            writer.writerow(info)
        print
        ptCount

        for filename in os.listdir(CaseFolder):
            ptID = filename.split('.')[0]
            ptFile = (os.path.join(CaseFolder, filename))
            if ptID in WTlist:
                status = 'WT'
                getInfo(ptFile, EA_WT, totalVarWT, EA_ALL, totalVarALL, geneListTotal)

            else:
                status = 'Case_Something_else'

            ptCount += 1
            print
            ptCount
            info = ([ptID, status])
            writer.writerow(info)

        print
        ptCount

    ##h5py
    sumEA_Risk = getSumEA(EA_Risk)
    sumEA_WT = getSumEA(EA_WT)
    sumEA_ALL = getSumEA(EA_ALL)

    save2h5py(totalVarRisk, 'totalVarRisk')
    save2h5py(totalVarWT, 'totalVarWT')
    save2h5py(totalVarALL, 'totalVarALL')
    save2h5py(sumEA_Risk, 'sumEA_Risk')
    save2h5py(sumEA_WT, 'sumEA_WT')
    save2h5py(sumEA_ALL, 'sumEA_ALL')

    print
    '-----Initializing RVEA calculation-----'
    # calculate RVEA of AA and GG using same fitted line (of all GG, AA, GA people)
    normRVEA_Risk, normRVEA_WT = calcRVEA_OneFit(geneListTotal, totalVarALL, sumEA_ALL, totalVarRisk, sumEA_Risk,
                                                 totalVarWT, sumEA_WT, '_real')
    save2h5py(normRVEA_Risk, 'RVEA_Risk')
    save2h5py(normRVEA_WT, 'RVEA_WT')

    RofRVEA = {}
    RofRVEA, RVEA_list = calcRofRVEA(geneListTotal, normRVEA_Risk, normRVEA_WT, RofRVEA, '_real')

    save2h5py(RofRVEA, 'RofRVEA')

    print
    '-----Initializing random RVEA calculation-----'
    ################
    #    RANDOM    #
    ################

    totalPt = Risklist + WTlist
    # germlineFiles = glob.glob(inputFolder + '/*.trauma')

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
            if i < len(Risklist):  ##Risk
                randAAlist.append(pt)
            else:
                randGGlist.append(pt)

        randAAlist_total.append(randAAlist)
        randGGlist_total.append(randGGlist)

    df = pd.DataFrame(randAAlist_total)
    df.to_csv(targetDirc + 'randomptsetsRisk.csv', sep='\t', index=False)

    df = pd.DataFrame(randGGlist_total)
    df.to_csv(targetDirc + 'randomptsetsWT.csv', sep='\t', index=False)

    # inputFolder = '/media/vision/ExtraDrive1/Exome/ALZ/ADSP/All_shortname/'
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
            if pt in Risklist:
                ptFile = (os.path.join(ControlFolder, pt))
            elif pt in WTlist:
                ptFile = (os.path.join(CaseFolder, pt))
            else:
                print
                'error'
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
                print
                '1000 random people done'

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

        print
        'random # = ', ptidx

    save2h5py(RVISofRVIS_random, 'RofRVEA_random')

    # Get z-scores
    noInfo = 0
    outputFile = targetDirc + 'ControlledZscores'
    with open(outputFile, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Gene', 'residual', 'mean', 'std', 'z-score', 'RofRVEA_random'])
        for gene in geneListTotal:
            try:
                mean = np.mean(RVISofRVIS_random[gene])
                std = np.std(RVISofRVIS_random[gene])
                x = RofRVEA[gene]
                z = (x - mean) / std
                r = RVISofRVIS_random[gene]
                new_r = []
                # for d in r:
                #    new_r.append(float('{0:.2f}'.format(d)))
                # info = ([gene, float('{0:.2f}'.format(x)), float('{0:.2f}'.format(mean)), float('{0:.2f}'.format(std)), float('{0:.2f}'.format(z)), new_r])

                info = ([gene, x, mean, std, z, r])
                writer.writerow(info)
            except KeyError:
                noInfo += 1

    print
    'number of genes with no z-score = ', noInfo


