'''
Created on Oct 18, 2019

@author: ywkim
'''

import scipy.stats as stats


def getOverlap(listR, listE, total):
    overlap = list(set(listR).intersection(listE))
    p = stats.hypergeom.sf(int(len(overlap)) - 1, int(total), int(len(listR)), int(len(listE)))

    return len(overlap), p, overlap


if __name__ == '__main__':
    #'''
    #overlap = 69
    overlap = 62
    listA = 166
    listB = 134
    total = 17039

    '''
    overlap = 22
    listA = 31
    listB = 31
    total = 7787
    '''

    p = stats.hypergeom.sf(int(overlap) - 1, int(total), int(listA), int(listB))

    print(p)