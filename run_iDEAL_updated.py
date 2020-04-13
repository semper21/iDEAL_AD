'''
Created on Jan 28, 2020
(Started pretty much re-writing on Apr 8, 2020)

@author: ywkim
'''

import pandas as pd
from pathlib import Path
from collections import Counter
from germline_analyses import get_matrix_as_df


def initial_lin_reg_one_fit():
    pass

if __name__ == '__main__':
    #input_folder = '/Users/ywkim/rosinante/ADSP/iDEAL_input_folder/' # This would change depending on
    input_folder = '/lab/rosinante/shared/ADSP/iDEAL_input_folder/'   # where you have Rosinante mounted on
    output_folder = str(Path().absolute()) + '/output/'  # For now the output files will be stored locally

    df_ADe2_sum = get_matrix_as_df(output_folder, 'sum_ea_matrix_ADe2', sep='\t') # TODO: needs to be changed to .tsv
    df_HCe4_sum = get_matrix_as_df(output_folder, 'sum_ea_matrix_HCe4', sep='\t')
    df_ADe2_freq = get_matrix_as_df(output_folder, 'frequency_matrix_ADe2', sep='\t')
    df_HCe4_freq = get_matrix_as_df(output_folder, 'frequency_matrix_HCe4', sep='\t')


