"""
Created on Aug 18, 2020

@author: ywkim
"""

import os
import csv
from pathlib import Path
from src.germline_analyses import get_list_from_csv


def get_variant_info(file_, gene_, variant_, flag):
    for line in open(file_):
        if line[0] == '#':
            continue
        cols = line.strip().split('\t')
        chr = str(int(float(cols[0])))
        pos = str(int(float(cols[1])))
        ref = cols[2]
        alt = cols[3]



        genes = cols[4].split(';')[0]
        sub = cols[5]

        if genes == gene_ and sub == variant_:
            sub_ = '/'.join((ref, alt))
            sift_ = [chr, pos, 1, sub_]

            chromosome = ''.join(('chr',chr))
            position = ':'.join((chromosome, pos))
            pph2_ = [position, sub_]

            cadd_ = [chr, pos, '.', ref, alt]

            if sift_ not in sift_list:
                sift_list.append(sift_)
                pph2_list.append(pph2_)
                cadd_list.append(cadd_)
                flag += 1
                break

    return flag


def output_variant_info(name, output_lists):
    output_file = output_folder + name + '_top_protective_variants.txt'
    with open(output_file, 'w') as f:
        if name == 'pph2':
            writer = csv.writer(f, delimiter=' ')
        else:
            writer = csv.writer(f, delimiter=',')

        for idx, output_list in enumerate(output_lists):
            gene = gene_list[idx]
            variant = variant_list[idx]
            # info = [gene, variant] + output_list
            info = output_list
            writer.writerow(info)


if __name__ == '__main__':
    input_folder = str(Path().absolute()) + '/input/'
    output_folder = str(Path().absolute()) + '/output_ADSP_discovery/'

    variant_file = input_folder + 'top_protective_variants.txt'
    # variant_file = input_folder + 'top_pathogenic_variants.txt'
    gene_list = get_list_from_csv(variant_file, 'Gene', sep='\t')
    variant_list = get_list_from_csv(variant_file, 'Sub', sep='\t')

    trauma_folder = '/media/vision/ExtraDrive1/Exome/ADSP_discovery/all_short_name/'

    sift_list = []
    pph2_list = []
    cadd_list = []
    for idx, gene in enumerate(gene_list):
        variant = variant_list[idx]
        flag = 0
        print(gene, variant)
        pt_counter = 0
        for filename in os.listdir(trauma_folder):
            ptFile = (os.path.join(trauma_folder, filename))
            flag = get_variant_info(ptFile, gene, variant, flag)
            pt_counter += 1
            if flag != 0:
                print(pt_counter)
                break



    # output_variant_info('sift', sift_list)
    # output_variant_info('pph2', pph2_list)
    output_variant_info('cadd', cadd_list)