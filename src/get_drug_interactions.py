'''
Created on May 7, 2020

@author: ywkim
'''

import pandas as pd
from pathlib import Path


if __name__ == '__main__':
    input_folder = str(Path().absolute()) + '/input/'
    output_folder = str(Path().absolute()) + '/output_ADSP_discovery/'

    ideal_file = input_folder + 'iDEAL_genelist.txt'
    df = pd.read_csv(ideal_file, sep='\t')

    dfDrug = pd.read_csv('/Users/ywkim/Desktop/Projects/GermlineProject/ADSP/RVEA/new2018/'
                         'RVEA_BaylorPass_nonHisWhite_2v4_h5py_STARTLOSS100/DGIdb_interactions_022119.tsv', sep='\t')
    dfDrug = dfDrug[['gene_name', 'drug_claim_primary_name', 'interaction_types', 'interaction_claim_source', 'PMIDs']]
    dfDrug.columns = ['Gene', 'Drug', 'InteractionType', 'Source', 'PMIDs']
    dfDrug = dfDrug.fillna('-')
    dfDG = df.merge(dfDrug, on=['Gene'], how='inner')

    dfDG.to_csv(output_folder + '216_drug_interactions.tsv', sep='\t', index=False)
