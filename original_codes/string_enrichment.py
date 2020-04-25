"""
Created on Oct 21, 2019

@author: ywkim
"""

import numpy as np
import pandas as pd
from random import sample
from sys import argv


def count_edges(test_list, truth_list, df_):
    """
    Count number of edges between two lists
    :param test_list: test gene list (either iDEAL hit genes or random 216 genes (negative controls))
    :param trust_list: gold standard gene list (AD genes)
    :param df: the network to use
    :return: an integer
    """

    df_subset = df_.loc[df_['protein1'].isin(test_list)]
    df_subset = df_subset.loc[df_subset['protein2'].isin(truth_list)]
    df_subset = df_subset.loc[df_subset['combined_score']>400]

    edges = len(set(df_subset['protein1'].values.tolist()))

    return edges


if __name__ == '__main__':
    #dirc = argv[1] # STRING_analysis/
    dirc = '/Users/ywkim/Desktop/Projects/GermlineProject/ADSP/RVEA/new2018/' \
           'RVEA_BaylorPass_nonHisWhite_2v4_h5py_STARTLOSS100/STRING_analysis/'

    """216 genes"""
    iDEAL_file = dirc.rsplit('/',2)[0] + '/iDEAL_genelist.txt'
    df_genes = pd.read_csv(iDEAL_file)
    hit_genes = df_genes['Gene'].values.tolist()

    """STRING v11.0"""
    network_file = dirc + '9606.protein.links.v11.0.txt'  # v11 downloaded on Oct 22, 2019 12:40PM
    df_network = pd.read_csv(network_file, sep=' ')

    """mapping ENSP to Gene symbol"""
    mapping_file = dirc + 'human.name_2_string.tsv'
    df_mapping = pd.read_csv(mapping_file, sep='\t', header=None, names=['taxid', 'Gene', 'ENSP'],
                             index_col=False)  # 19095

    map_dict = dict(zip(df_mapping.ENSP, df_mapping.Gene))  # 19095

    df_network['protein1'] = df_network['protein1'].map(map_dict)
    df_network['protein2'] = df_network['protein2'].map(map_dict)

    all_gene_file = dirc.rsplit('/', 2)[0] + '/ControlledZscores'
    df_all_genes = pd.read_csv(all_gene_file, sep='\t')
    all_genes = df_all_genes['Gene'].values.tolist()

    n = len(hit_genes)
    df = pd.DataFrame(columns=[i for i in range(1000)])
    for i in range(1000):
        random_draw = sample(all_genes, int(n))
        df[i] = random_draw

    df.to_csv(dirc + '1000x_random_genes_to_get_edges.csv', sep=',', index=False)

    """GWAS genes"""
    gene_types = ['BaitGWAS', 'BaitGWAS_FAD', 'APP_Tau']
    truth_file = dirc.rsplit('/',2)[0] + '/216_updated/genes_for_network_anaylsis.xlsx'
    for gene_type in gene_types:
        df_gwas = pd.read_excel(truth_file, sheet_name=gene_type)
        gwas_genes = df_gwas['Gene'].values.tolist()

        edge_number = count_edges(hit_genes, gwas_genes, df_network)

        random_edges_list = []
        for i in range(1000):
            if i % 100 == 0:
                print(i)
            random_genes = df[i].values.tolist()
            random_edge_number = count_edges(random_genes, gwas_genes, df_network)
            random_edges_list.append(random_edge_number)

        df_edges = pd.DataFrame(random_edges_list, columns = ['RandomEdges'])
        df_edges.to_csv(dirc + gene_type + '_random_edges_distribution.tsv', sep='\t', index=False)

        mean = np.mean(random_edges_list)
        std = np.std(random_edges_list)

        z = (edge_number - mean) / std

        print(gene_type, edge_number, z)

