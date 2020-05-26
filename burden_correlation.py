'''
Created on May 12, 2020

@author: ywkim
'''
import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt

from germline_analyses import get_matrix_subset, get_list_from_csv

from sklearn.cluster import DBSCAN
from sklearn.cluster import SpectralClustering
from sklearn.cluster import KMeans

if __name__ == '__main__':
    input_folder = str(Path().absolute()) + '/input/'
    output_folder =  str(Path().absolute()) + '/output/'

    phenotype_file = input_folder + 'ADSP_phenotype.tsv'

    gene_file = output_folder + 'iDEAL_genelist.txt'
    gene_list = get_list_from_csv(gene_file, 'Gene', sep='\t')

    ADe2 = get_matrix_subset(output_folder, 'sum_ea_matrix_ADe2.tsv', gene_list, sep='\t', index_col=0)
    ADe2_sum = ADe2.sum(axis=0)

    patient_list = pd.read_csv(output_folder + 'sum_ea_matrix_ADe2.tsv', sep='\t', index_col=0).columns.tolist()

    df = pd.read_csv(phenotype_file, sep='\t', index_col=None)
    age_list = []
    for patient in patient_list:
        age = df.loc[df['SUBJID'] == patient, 'Age'].values[0]
        age_list.append(age)

    df_survival = pd.DataFrame({'Patient': patient_list, 'sum_EA':ADe2_sum, 'Age': age_list})
    df_survival.to_csv(output_folder + 'survival_sumEA_age.tsv', sep='\t', index=False)
