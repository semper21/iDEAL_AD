'''
Created on May 5, 2020

@author: ywkim
'''

from sys import argv
from germline_analyses import get_list_from_csv, enrichment_test


if __name__ == '__main__':

    input_file1 = argv[1]
    input_file2 = argv[2]

    l1 = get_list_from_csv(input_file1, 'Gene', sep='\t')
    l2 = get_list_from_csv(input_file2, 'Gene', sep='\t')

    # n_total = len(list(set(l1+l2)))
    n_overlap, p, overlap = enrichment_test(l1, l2, 17039)

    print(len(l1), len(l2))
    print(n_overlap, p)
    print(overlap)

