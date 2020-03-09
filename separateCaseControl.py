'''
Created on Jul 25, 2019

@author: ywkim
'''
import os
import glob

def processPhenoFile(input_file):
    excluded = 0
    caseList = []
    controlList = []
    for line in open(input_file):
        if line[0] == '#' or line[0] == 'd':
            continue
        cols = line.strip().split('\t')
        # print cols
        subID = cols[1]
        state = cols[13]

        if state == '0' or state == 0: #Healthy controls
            controlList.append(subID)
        elif state == '1' or state == 1: #AD patients
            caseList.append(subID)
        else:
            excluded += 1 #Most likely people misdiagnosed

    print(excluded)
    return caseList, controlList

if __name__ == '__main__':
    target_directory = '/lab/rosinante/shared/ADSP/iDEAL_input_folder/' #directory on rosinante
    #target_directory = '/media/vision/ExtraDrive1/Exome/ADSP/' #my local machine
    phenotype_file = target_directory + 'phs000572.v7.pht005179.v1.p4.c1.CaseControlEnrichedPhenotypesWES_y1.HMB-IRB.txt'

    caseList, controlList = processPhenoFile(phenotype_file)

    srcglob = glob.glob(target_directory + 'after_QCfilter_jamie/*.trauma')

    os.mkdir(target_directory + 'Case/')
    os.mkdir(target_directory + 'Control/')
    os.mkdir(target_directory + 'all_short_name/')  # This was needed to more conveniently match patient IDs from
                                                    # the phenotype file

    for f in range(len(srcglob)):
        src = srcglob[f]  # each file
        pt = src.rsplit('/', 1)[1].split('-')[0:3]
        pt = '-'.join(pt)
        print(pt)
        if pt in caseList:
            os.symlink(src, target_directory + 'Case/' + pt)
        elif pt in controlList:
            os.symlink(src, target_directory + 'Control/' + pt)
        else:
            print('Patient neither case nor control')

        os.symlink(src, target_directory + 'all_short_name/' + pt)
    print(str(len(caseList) + len(controlList)) + ' out of 5686')
