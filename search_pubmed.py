'''
Created on Aug 14, 2020

@author: ywkim
'''

import csv
import time
import pandas as pd
from Bio import Entrez
from pathlib import Path
from germline_analyses import get_list_from_csv
from urllib.error import HTTPError


def search(query, rs):
    term_ = 'Alzheimer AND ' + query
    try:
        Entrez.email = 'ywkim@bcm.edu'
        handle = Entrez.esearch(db='pubmed',
                                sort='relevance',
                                retstart=rs,
                                retmax='100000',
                                retmode='xml',
                                term=term_)
        results = Entrez.read(handle)

    except HTTPError:
        time.sleep(5)
        Entrez.email = 'ywkim@bcm.edu'
        handle = Entrez.esearch(db='pubmed',
                                sort='relevance',
                                retstart=rs,
                                retmax='100000',
                                retmode='xml',
                                term=term_)
        results = Entrez.read(handle)

    return results

if __name__ == '__main__':
    input_folder = str(Path().absolute()) + '/input/'
    output_folder = str(Path().absolute()) + '/output/'

    drugs_file = input_folder + 'iDEAL_drugs_list.txt'
    drugs_list = get_list_from_csv(drugs_file, 'Drug', sep='\t')

    outputFile = output_folder + 'iDEAL_drugs_pubmed.txt'
    with open(outputFile, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Drug', 'PMID'])
        for idx, drug in enumerate(drugs_list):
            print(drug)
            id_list = search(drug, 0)['IdList']

            writer.writerow([drug, id_list])
