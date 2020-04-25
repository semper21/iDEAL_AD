"""
Created on March 9, 2020

@author: ywkim
"""

import os
import pandas as pd
from pathlib import Path

if __name__ == '__main__':
    # input_folder = '/Users/ywkim/rosinante/shared/ADSP/iDEAL_input_folder/' # this would change depending on
    # input_folder = '/lab/rosinante/shared/ADSP/iDEAL_input_folder/'  # where you have Rosinante mounted on
    input_folder = '/lab/rosinante/shared/ADSP/iDEAL_input_folder/ADSP_extension/'
    output_folder = str(Path().absolute()) + '/output_ADSP_extension/'  # for now the output files will be stored locally

    """
    Need to get the list of all the genes sequenced/observed
    """
    total_gene_set = set()
    #trauma_folder = input_folder + 'after_QCfilter_jamie/'
    trauma_folder = input_folder + 'Actions/'

    for idx, filename in enumerate(os.listdir(trauma_folder)):
        pt_file = os.path.join(trauma_folder, filename)
        df = pd.read_csv(pt_file, sep='\t', index_col=None)
        df.GENE = df.GENE.astype(str)
        genes = list(set(df['GENE'].values.tolist()))
        new_genes = [gene.split(';')[0] for gene in genes if gene != 'nan']

        total_gene_set.update(new_genes)
        print(idx)

    total_gene_list = list(sorted(total_gene_set))
    print(len(total_gene_list))

    df_genes = pd.DataFrame({'Gene':total_gene_list})
    df_genes.to_csv(output_folder + 'total_gene_list.csv', sep=',', index=False)