"""
Created on 

@author: ywkim

#

"""
import os
import csv
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from sklearn import metrics
from sklearn.metrics import confusion_matrix
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split



def get_quality(file_):
    dict_ = {}
    for line in open(file_):
        if line[0] == '#':
            continue
        cols = line.strip().split('\t')
        chs = cols[0]
        pos = cols[1]
        ref = cols[2]
        alt = cols[3]

        filt = cols[4]

        s = '-'
        mut = (chs, pos, ref, alt)
        substitution = s.join(mut)

        dict_[substitution] = filt

    return dict_


def get_phenotype(phenotype_file_):
    protective = []
    risk = []
    neutral_hc = []
    neutral_ad = []
    for line in open(phenotype_file_):
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
                    risk.append(sub_id)     # healthy with risk
                if apoe == '33':
                    neutral_hc.append(sub_id)

            elif state == '1' or state == 1:
                if apoe == '22' or apoe == '23':    # AD with protection
                    protective.append(sub_id)
                if apoe == '33':
                    neutral_ad.append(sub_id)
            else:
                pass

    return protective, risk, neutral_hc, neutral_ad


def convert_strings(ea_list):
    """ For variants affecting multiple isoforms and therefore have multiple EA annotations,
        converts any string to float

    :param ea_list: list of EA scores of all the isoforms affected by the variant
    :return: a (temporary) list of EA scores without any strings
    """
    temp_list = []
    for x in ea_list:
        try:
            temp_list.append(float(x))
        except ValueError:  # When there is a string
            if x in ['STOP', 'no_STOP', 'STOP-loss', 'START_loss']:
                temp_list.append(100.0)
            else:
                continue
    return temp_list


def process_germline_file(trauma_file_):
    ea_dict = {}
    for line in open(trauma_file_):
        if line[0] == '#':
            continue
        cols = line.strip().split('\t')
        chs = cols[0]
        pos = cols[1]
        ref = cols[2]
        alt = cols[3]
        gene = cols[4].split(';')[0]
        sub = cols[5]
        action = cols[6]

        if gene == '':
            continue

        s = '-'
        mut = (chs, pos, ref, alt)
        substitution_ = s.join(mut)

        try:
            quality = qc_dict[substitution_]
        except KeyError:
            quality = 'non-PASS'

        if quality == 'non-PASS':
            continue
        elif quality == 'PASS':
            pass
        else:
            print('ERROR')

        if sub == '-':      # no indels
            continue
        if gene == 'APOE' and sub == 'C130R':
            continue

        if gene in gene_list:
            if gene not in ea_dict:
                ea_dict[gene] = []
            if action in ['silent', '-', 'no_action', 'no_trace', 'no_gene'] or gene == '':
                continue
            elif action in ['STOP', 'no_STOP', 'STOP-loss', 'START_loss']:
                ea_dict[gene].append(100.)
            else:
                try:
                    spliced_ea = action.strip().split(';')
                    new_spliced_ea = convert_strings(spliced_ea)
                    if len(new_spliced_ea) == 0:
                        pass
                    else:
                        average = np.mean(new_spliced_ea)
                        ea_dict[gene].append(average)
                except AttributeError:
                    ea_dict[gene].append(float(action))

    return ea_dict


def fill_matrix(trauma_file_, matrix_, pt_idx):
    ea_dict = process_germline_file(trauma_file_)
    for gene_idx, gene_ in enumerate(gene_list):
        try:
            ea_list = ea_dict[gene_]
            if len(ea_list) != 0:
                sum_ea = np.sum(ea_list)
                matrix_[gene_idx][pt_idx] = sum_ea / len(ea_list)
            else:
                matrix_[gene_idx][pt_idx] = 0

        except KeyError:
            matrix_[gene_idx][pt_idx] = 0


def convert_to_matrix(exome_folder_, status_list):
    pt_count = 0
    pt_idx = 0
    matrix = np.zeros((len(gene_list), len(status_list)))
    for filename in os.listdir(exome_folder_):
        pt_id = filename.split('.')[0]
        pt_file = (os.path.join(exome_folder_, filename))
        if pt_id in status_list:
            fill_matrix(pt_file, matrix, pt_idx)
            pt_idx += 1
        else:
            pass
        pt_count += 1
        print(pt_count)

    return matrix


if __name__ == '__main__':
    # target_directory = '/home/vision/Documents/GermlineProject/ADSP/RVEA_BaylorPass_nonHisWhite_2v4_h5py_STARTLOSS100/'
    target_directory = '/home/vision/Documents/GermlineProject/ADSP/iDEAL_updated_test/'
    exome_folder = '/media/vision/ExtraDrive1/Exome/ADSP/'

    gene_file = target_directory + 'iDEAL_genelist.txt'
    df_gene = pd.read_csv(gene_file, index_col=None)
    gene_list = df_gene['Gene'].values.tolist()

    quality_file = '/home/vision/Documents/GermlineProject/ADSP/snvquality_detailed_jamie.csv'
    phenotype_file = exome_folder + 'phs000572.v7.pht005179.v1.p4.c1.CaseControlEnrichedPhenotypesWES_y1.HMB-IRB.txt'

    case_folder = exome_folder + 'Case/'
    control_folder = exome_folder + 'Control/'

    qc_dict = get_quality(quality_file)

    protective_list, risk_list, neutral_hc_list, neutral_ad_list = get_phenotype(phenotype_file)

    matrix_AD_protective = convert_to_matrix(case_folder, protective_list)
    matrix_HC_risk = convert_to_matrix(control_folder, risk_list)
    matrix_AD_neutral = convert_to_matrix(case_folder, neutral_ad_list)
    matrix_HC_neutral = convert_to_matrix(control_folder, neutral_hc_list)


    """Logistic Regression"""
    x = np.vstack((matrix_AD_protective.T, matrix_HC_risk.T))
    y = [1] * len(matrix_AD_protective.T) + [0] * len(matrix_HC_risk.T)

    x_train, x_test, y_train, y_test = train_test_split(x, y, test_size = 0.2)

    lr = LogisticRegression(random_state=1)
    lr.fit(x_train, y_train)

    y_pred = lr.predict(x_test)
    cfm = confusion_matrix(y_test, y_pred)
    print('test confusion matrix = ', cfm)

    print("Accuracy:", metrics.accuracy_score(y_test, y_pred))
    print("Precision:", metrics.precision_score(y_test, y_pred))
    print("Recall:", metrics.recall_score(y_test, y_pred))

    x_val = np.vstack((matrix_AD_neutral.T, matrix_HC_neutral.T))
    y_val = [1] * len(matrix_AD_neutral.T) + [0] * len(matrix_HC_neutral.T)

    y_pred_val = lr.predict(x_val)
    cfm2 = confusion_matrix(y_val, y_pred_val)
    print('validation confusion matrix = ', cfm2)

    plt.clf()
    sns.heatmap(cfm2, annot=True, cmap='Spectral_r')

    print("Accuracy:", metrics.accuracy_score(y_val, y_pred_val))
    print("Precision:", metrics.precision_score(y_val, y_pred_val))
    print("Recall:", metrics.recall_score(y_val, y_pred_val))

    plt.show()