'''
Created on Apr 22, 2020

@author: ywkim
'''

import zarr
import allel
import numpy as np
import pandas as pd

#TODO: THIS NEEDS TO BE RE-RUN
def get_phenotype_from_excel(phenotype_file, samples_):
    df_= pd.read_excel(phenotype_file)

    samples_callset_index = [samples_.index(s) for s in df_['ID']]
    df_['callset_idx'] = samples_callset_index

    all_pt = df_['ID'].values.tolist()
    all_idx = df_['callset_idx'].values.tolist()
    HC2 = df_.loc[(df_['AD'] == 0) & (df_['APOE'].isin([22, 23]))]['ID'].values.tolist()
    HC3 = df_.loc[(df_['AD'] == 0) & (df_['APOE'] == 33)]['ID'].values.tolist()
    HC4 = df_.loc[(df_['AD'] == 0) & (df_['APOE'].isin([44, 34]))]['ID'].values.tolist()
    AD2 = df_.loc[(df_['AD'] == 1) & (df_['APOE'].isin([22, 23]))]['ID'].values.tolist()
    AD3 = df_.loc[(df_['AD'] == 1) & (df_['APOE'] == 33)]['ID'].values.tolist()
    AD4 = df_.loc[(df_['AD'] == 1) & (df_['APOE'].isin([44, 34]))]['ID'].values.tolist()

    HC2_idx = df_.loc[(df_['AD'] == 0) & (df_['APOE'].isin([22, 23]))]['callset_idx'].values.tolist()
    HC3_idx = df_.loc[(df_['AD'] == 0) & (df_['APOE'] == 33)]['callset_idx'].values.tolist()
    HC4_idx = df_.loc[(df_['AD'] == 0) & (df_['APOE'].isin([44, 34]))]['callset_idx'].values.tolist()
    AD2_idx = df_.loc[(df_['AD'] == 1) & (df_['APOE'].isin([22, 23]))]['callset_idx'].values.tolist()
    AD3_idx = df_.loc[(df_['AD'] == 1) & (df_['APOE'] == 33)]['callset_idx'].values.tolist()
    AD4_idx = df_.loc[(df_['AD'] == 1) & (df_['APOE'].isin([44, 34]))]['callset_idx'].values.tolist()

    return all_pt, all_idx, AD2, AD3, AD4, HC2, HC3, HC4, HC2_idx, HC3_idx, HC4_idx, AD2_idx, AD3_idx, AD4_idx

if __name__ == '__main__':
    dirc = '/Users/ywkim/rosinante/ADSP/iDEAL_input_folder/ADSP_extension/'

    phenotype_file = dirc + 'ADSP_extension_phenotypes.xlsx'
    vcf_file = dirc + '/anno.coding.4.gcad.qc.wgs.4789.GATK.2018.09.17.biallelic.genotypes.ALL.snps.liftover.chr.' \
                      'sorted.pass.reheader3.noCHR.whites.nonzero.vcf.bz2.EA.VCF'   # TODO: this needs to be changed

    zarr_file = vcf_file + '.zarr'
    allel.vcf_to_zarr(vcf_file, zarr_file, fields=['samples', 'calldata/GT', 'variants/gene', 'variants/sub', 'variants/EA'])

    callset = zarr.open_group(zarr_file, mode='r')

    samples = list(callset['samples'][:])

    total_pt_list, total_pt_idx, ADe2, ADe3, ADe4, HCe2, HCe3, HCe4, \
    ADe2_idx, ADe3_idx, ADe4_idx, HCe2_idx, HCe3_idx, HCe4_idx = get_phenotype_from_excel(phenotype_file, samples)

    ea_list = list(callset['variants/EA'][:])
    sub_list = list(callset['variants/gene'][:])
    gene_list = list(callset['variants/gene'][:])

    total_gene_list = list(set(gene_list))  # TODO: This will have to be output-ed

    gt_zarr = callset['calldata/GT']
    gt = allel.GenotypeArray(gt_zarr)   # (n_variant, n_sample, ploidy)
    # gt_matrix = np.asarray(gt.to_packed())  # 0 = 0/0, 1 = 0/1, 17 = 1/1, 239 = -1/-1

    pt_lists = [ADe2, ADe3, ADe4, HCe2, HCe3, HCe4]
    pt_idx_lists = [ADe2_idx, ADe3_idx, ADe4_idx, HCe2_idx, HCe3_idx, HCe4_idx]
    for idx, group in enumerate(['ADe2', 'ADe3', 'ADe4', 'HCe2', 'HCe3', 'HCe4']):
        pt_list = pt_lists[idx]
        pt_idx = pt_idx_lists[idx]
        matrix_freq = np.zeros((len(total_gene_list), len(pt_list)))
        matrix_sum = np.zeros((len(total_gene_list), len(pt_list)))

        gt_subgroup = gt.take(pt_idx, axis=1)
        gt_sub_matrix = np.asarray(gt.to_packed())


