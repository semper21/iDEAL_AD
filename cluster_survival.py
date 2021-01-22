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
    output_folder =  str(Path().absolute()) + '/output_ADSP_discovery/'

    phenotype_file = input_folder + 'ADSP_phenotype.tsv'

    gene_file = output_folder + 'iDEAL_genelist.txt'
    gene_list = get_list_from_csv(gene_file, 'Gene', sep='\t')

    ADe2 = get_matrix_subset(output_folder, 'sum_ea_matrix_ADe2.tsv', gene_list, sep='\t', index_col=0)
    ADe2_sum = ADe2.sum(axis=0)
    patient_list = pd.read_csv(output_folder + 'sum_ea_matrix_ADe2.tsv', sep='\t', index_col=0).columns.tolist()

    # clustering = DBSCAN(eps=3, min_samples=10).fit(ADe2.T)
    # clustering = SpectralClustering(n_clusters=2, assign_labels="discretize", random_state=0).fit_predict(ADe2)

    clustering = KMeans(n_clusters=2).fit_predict(ADe2.T)

    idx_list1 = [i for i, x in enumerate(clustering) if x == 1]
    idx_list0 = [i for i, x in enumerate(clustering) if x == 0]

    pt_list1 = [patient_list[i] for i in idx_list1]
    pt_list0 = [patient_list[i] for i in idx_list0]

    df = pd.read_csv(phenotype_file, sep='\t', index_col=None)

    df1 = df[df.SUBJID.isin(pt_list1)]
    df0 = df[df.SUBJID.isin(pt_list0)]

    age1 = df1['Age'].values.tolist()
    age0 = df0['Age'].values.tolist()
    age10 = age1 + age0
    group = [1] * len(age1) + [0] * len(age0)

    df_survival = pd.DataFrame({'Group': group, 'Age': age10})
    df_survival.to_csv(output_folder + 'age_survival_kmeans_n_2.tsv', sep='\t', index=False)
