# Imputed Deviation of Evolutionary Action Load (iDEAL) analysis of Alzheimer's disease


## Project Description
Despite considerable insight from GWAS and sequencing studies, much of late-onset Alzheimer’s Disease (LOAD) heritability remains unexplained leading to uncertain risk stratification of patients and no available disease-modifying therapies. To date, the APOE gene remains the strongest genetic risk factor for LOAD. Carriers of the APOEɛ4 allele are at a greater risk, while the APOEɛ2 allele plays a protective role. However, many APOEɛ4 carriers remain disease-free and some APOEɛ2 carriers develop LOAD.

iDEAL, or imputed Deviation of Evolutionary Action Load, is a novel approach that quantifies the differential mutational burden for each gene between the two paradoxical patient groups through a series of linear regressions. 

The genes identified by iDEAL were significantly enriched for differentially expressed genes (DEGs) in the AD vs healthy control brains and were highly connected to known LOAD genes identified by GWAS. Furthermore, many of the genes were relevant in vivo as they ameliorated neurodegeneration caused by tau and secreted β42 using well-validated Drosophila models. Network analyses revealed involvement of candidate genes in brain cell-type specific pathways including synaptic biology, dendritic spine pruning and inflammation. Finally, in a 5-fold cross-validation, the AdaBoost-SVM algorithm trained on the mutational features across the iDEAL genes could predict which of the APOEɛ2 carriers would develop AD, and conversely, which APOEɛ4 carriers would remain healthy with average AUCs of 0.79 and 0.71, respectively.
## Publication

[Kim YW, Al-Ramahi I, Koire A, Wilson SJ, Konecki DM, Mota S, Soleimani S, Botas J, Lichtarge O.
Harnessing the Paradoxical Phenotypes of APOE2 and APOE4 to Identify Genetic Modiers in
Alzheimer's Disease. Alzheimer's & Dementia. 2020](https://alz-journals.onlinelibrary.wiley.com/doi/10.1002/alz.12240)


---
## Installation

### Download code
```
git clone https://github.com/semper21/iDEAL_AD.git 
```
### Install environment

Requirement = python 3.6

```
conda env create -f environment.yml
```