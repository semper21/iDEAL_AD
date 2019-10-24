"""
Created on Oct 24, 2019

@author: ywkim

Ismael wants to only use GWAS genes that interact with IDEAL genes as "positive control."
"""

import numpy as np
import pandas as pd
from random import sample
from sys import argv

from IPython import embed


def count_edges(test_list, truth_list, df_, cutoff, gene_type):
    """
    Count number of edges between two lists
    :param test_list: test gene list (either iDEAL hit genes or random 216 genes (negative controls))
    :param trust_list: gold standard gene list (AD genes)
    :param df: the network to use
    :param cutoff: minimum required interaction score
    :param gene_type: the type of trust list
    :return: # of genes in test_list that interact with genes in truth list (integer)

    The equation for combining scores of each channel if from: http://version11.string-db.org/help/faq/
    """

    df_subset = df_.loc[df_['protein1'].isin(test_list)]
    df_subset = df_subset.loc[df_subset['protein2'].isin(truth_list)]

    if gene_type == 'Bait3':
        df_subset = df_subset.loc[df_subset['combined_score'] > cutoff]
    else:
        p = 0.041
        '''If score < prior, we will end up with negative combined score - so we take care of that first'''
        channels = ['textmining', 'experimental', 'database']
        df_subset[channels] = np.where(df_subset[channels] < (p*1000), (p*1000), df_subset[channels])
        df_subset['txt'] = ((df_subset['textmining']*0.001) - p) / (1 - p)
        df_subset['exp'] = ((df_subset['experimental']*0.001) - p) / (1 - p)
        df_subset['dtb'] = ((df_subset['database']*0.001) - p) / (1 - p)
        df_subset['tot'] = 1 - ((1 - df_subset['txt']) * (1 - df_subset['exp']) * (1 - df_subset['dtb']))
        df_subset['final_score'] = (df_subset['tot'] + p * (1 - df_subset['tot'])) * 1000
        # embed()
        df_subset = df_subset.loc[df_subset['final_score'] > cutoff]

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
    network_file = dirc + '9606.protein.links.detailed.v11.0.txt'  # downloaded on Oct 24, 2019 12:10PM
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

    df.to_csv(dirc + '1000x_random_genes_to_get_edges2.csv', sep=',', index=False)


    """GWAS genes"""
    gene_types = ['Bait1', 'Bait2', 'Bait3', 'Bait4']
    confidences = [400, 450, 700, 400]
    truth_file = dirc.rsplit('/',2)[0] + '/216_updated/Bait_GWAS_genes.xlsx'
    for idx, gene_type in enumerate(gene_types):
        df_gwas = pd.read_excel(truth_file, sheet_name=gene_type)
        gwas_genes = df_gwas['Gene'].values.tolist()
        confidence = confidences[idx]
        edge_number = count_edges(hit_genes, gwas_genes, df_network, confidence, gene_type)

        random_edges_list = []
        for i in range(1000):
            if i % 100 == 0:
                print(i)
            random_genes = df[i].values.tolist()
            random_edge_number = count_edges(random_genes, gwas_genes, df_network, confidence, 'random')
            random_edges_list.append(random_edge_number)

        df_edges = pd.DataFrame(random_edges_list, columns = ['RandomEdges'])
        df_edges.to_csv(dirc + gene_type + '_random_edges_distribution.tsv', sep='\t', index=False)

        mean = np.mean(random_edges_list)
        std = np.std(random_edges_list)

        z = (edge_number - mean) / std

        print(gene_type, edge_number, z)

