'''
Created on Apr 21, 2020

@author: ywkim
'''

import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles
from src.germline_analyses import get_list_from_csv

def venn_diagram(list1, list2, label1, label2, title, filename):
    venn = venn2([set(list1), set(list2)], set_labels=(label1, label2), set_colors=('red', 'skyblue'), alpha=0.7)
    venn2_circles([set(list1), set(list2)])
    for text in venn.set_labels:
        text.set_fontsize(14)
    for text in venn.subset_labels:
        text.set_fontsize(16)
    plt.title(title)
    plt.savefig(new_ideal_folder + filename + '.png', dpi=200, transparent=True)
    plt.clf()


if __name__ == '__main__':
    new_ideal_folder = '/Users/ywkim/Desktop/Projects/GermlineProject/ADSP/updated_iDEAL/'
    old_ideal_folder = '/Users/ywkim/Desktop/Projects/GermlineProject/ADSP/RVEA/new2018/RVEA_BaylorPass_' \
                       'nonHisWhite_2v4_h5py_STARTLOSS100/'

    new_pathogenic_list = get_list_from_csv(new_ideal_folder + 'pathogenic.tsv', 'Gene', sep='\t')
    new_protective_list = get_list_from_csv(new_ideal_folder + 'protective.tsv', 'Gene', sep='\t')
    old_pathogenic_list = get_list_from_csv(old_ideal_folder + 'pathogenic.txt', 'Gene', sep='\t')
    old_protective_list = get_list_from_csv(old_ideal_folder + 'protective.txt', 'Gene', sep='\t')

    #'''
    venn_diagram(old_pathogenic_list, new_pathogenic_list, 'Old ver', 'New ver', 'APOE2-AD genes',
                      'APOE2_AD_genes_comparison')
    venn_diagram(old_protective_list, new_protective_list, 'Old ver', 'New ver', 'APOE4-HC genes',
                      'APOE4_HC_genes_comparison')
    #'''
