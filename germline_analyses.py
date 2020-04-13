'''
Created on Apr 6, 2020

@author: ywkim
'''

import pandas as pd

def delNoStr(current_list):
    temp_list = []
    for ea in current_list:
        try:
            temp_list.append(float(ea))
        except ValueError:  # When there is a string
            if ea in ['STOP', 'no_STOP', 'STOP-loss', 'START_loss']:
                temp_list.append(1.0)
            else:
                continue
    return temp_list


def matrix_out(matrix_, pt_list_, gene_list_, output_dirc, name):
    """ Outputs matrices to files

    :param matrix_: matrix
    :param pt_list_: patient list
    :param gene_list_: gene list
    :param output_dirc: output directory
    :param name: file name
    :return: -
    """
    df = pd.DataFrame(matrix_, index=gene_list_)
    df.columns = pt_list_
    df.to_csv(output_dirc + name + '.tsv', sep='\t', index=True)


def output_dict(any_dict, output_dirc, name, sep):
    """ Outputs any dictionary to a file

    :param any_dict: any dictionary
    :param output_dirc: output directory
    :param name: file name
    :return: -
    """
    df = pd.DataFrame.from_dict(any_dict, orient='index')
    df.to_csv(output_dirc + name, sep=sep, index=True)


def get_matrix_as_df(matrix_folder, filename, sep):
    matrix_file = matrix_folder + filename
    df = pd.read_csv(matrix_file, sep=sep, index_col=None)

    return df


if __name__ == '__main__':
    pass