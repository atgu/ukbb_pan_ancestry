---
id: qc
title: Quality Control (QC)
---

## Phenotypes

- We ran PHESANT on continuous (including biomarker) and categorical phenotypes. This performed inverse-rank normal transformation (IRNT) of all continuous traits within each ancestry group.
- We did not run PHESANT on ICD-10 codes or prescription data.
- We constructed [Phecodes](https://phewascatalog.org/phecodes_icd10) from ICD-9 and ICD-10 codes using [our custom script](https://github.com/atgu/ukbb_pan_ancestry/blob/master/assign_phecodes.py) based on the interim outputs of the [createUKBphenome](https://github.com/umich-cphds/createUKBphenome) scripts. ICD-9 codes include primary codes (data field: 41203) and secondary codes (41205), and ICD-10 codes include primary codes (41202), secondary codes (41204), external codes (41201), and cause of death codes (40001). The mappings between phecodes and ICD codes are summarized [here](https://github.com/atgu/ukbb_pan_ancestry/blob/master/data/UKB_Phecode_v1.2b1_ICD_Mapping.txt).
- For binary traits, we required at least 50 cases in each population except EUR given the larger sample size. In EUR, we required at least 100 cases.
- Prescription data was extracted from UKB data field 42039 and processed through a custom pipeline to harmonize data.
- Special treatment of some phenotypes was applied, including for those related to COVID-19 as well as custom phenotype combinations such as waist-hip ratio (phenotype computed from UKB code 48 / 49, then IRNT).

Pre-QC phenotype counts are broken down by ancestries as follows:

| Population   | Total phenotypes | Categorical | Continuous | Phecode | ICD-10 | Biomarkers | Prescriptions |
|-------|----------------|------------|---------|-------|-------------|------------|---------------|
| AFR |          2493 |        981 |     337 |   197 |        725  |         30 |           223 |
| AMR |          1105 |        423 |      31 |    20 |        561  |         30 |            40 |
| CSA |          2771 |       1051 |     418 |   234 |        719  |         30 |           319 |
| EAS |          1612 |        618 |      91 |    55 |        714  |         29 |           105 |
| EUR |          7200 |       3672 |    1325 |   929 |        800  |         30 |           444 |
| MID |          1372 |        509 |      83 |    52 |        591  |         30 |           107 |


## Sample QC

- We removed participants flagged with `putative.sex.chromosome.aneuploidy == 1`

### Ancestry definitions

We first combined reference data from the 1000 Genomes Project and Human Genome Diversity Panel (HGDP). Briefly, we used 1000 Genomes phase 3 sequenced data and publicly available genotyping data from HGDP. We combined these reference datasets into continental ancestries according to their corresponding meta-data. Details of meta-data grouping combining these two datasets are available [here](https://docs.google.com/spreadsheets/d/1jenSz_HnbA1kBESaUmur3Ob72-EPXJgfUWhbz5UdltA/edit#gid=433808438).

We then used a two-stage approach to: 
1. assign continental ancestries, and then
2. prune ancestry outliers within continental groups. 

First, we ran PCA on unrelated individuals from the 1000 Genomes + HGDP combined reference dataset. To partition individuals in the UK Biobank based on their continental ancestry, we used the PC loadings from the reference dataset to project UK Biobank individuals into the same PC space. We trained a random forest classifier given continental ancestry meta-data based on top 6 PCs from the reference training data. (We also assessed whether including more PCs -- up to 20 -- changed the ancestry assignments but identified no change, so we continued with the assignments from the top 6 PCs). We applied this random forest to the projected UK Biobank PCA data and assigned initial ancestries that were subsequently refined if the random forest probability was >50%. Other individuals with a probability < 50% for any given ancestry group were dropped from further analysis. 

Second, we refined initial ancestry assignments by pruning outliers within each continental assignment. We started by rerunning PCA among UK Biobank individuals within each assigned continental ancestry group (i.e. excluding reference panel data). We calculated the total distance from population centroids across 10 PCs. Using the PC scores, we computed centroid distances across 3-5 centroids spanning these PCs depending on the degree of heterogeneity within each continental ancestry as follows:

For each individual and for each ellipse across 10 PCs, subtract the population PC mean from the individual's PC, square this value, and divide by the variance of the PC.

We identified ancestry outliers based by plotting histograms of centroid distances and removing those individuals from the extreme high end of the distribution.

### Relatedness

Within each population, we ran PC-Relate with `k=10` and `min_individual_maf=0.05`. To get the maximal set of unrelateds, we then ran `hl.maximal_independent_set` using Hail. We used these unrelated individuals to define the ancestry-specific PC space, then projected related individuals into the same PC space for use as covariates in SAIGE.

## Variant QC

- We used imputed variants from the UK Biobank (97,059,328 in version 3). 
- We retained variants with INFO scores > 0.8 (29,865,259 variants on the autosomes and X chromosome). - For each population, we retained variants with an allele count of at least 20, computing using the genotypes where available or imputed dosages.
  - We created a low-confidence filter for those variants with a <b><u>minor</u></b> allele count less than 20.

## GWAS model

We ran GWAS for all phenotypes and ancestry groups using the Scalable and Accurate Implementation of GEneralized mixed model package, [SAIGE](https://www.nature.com/articles/s41588-018-0184-y), which runs a linear or logistic mixed model including a kinship matrix as a random effect and [covariates](https://github.com/atgu/ukbb_pan_ancestry/wiki/QC/_edit#covariates) as fixed effects. This approach was chosen due to its suitability for a wide range of GWAS scenarios, its ability to account for more heterogeneity in sample makeup, and its computational scalability. Specifically, it provides accurate p values even when case-control ratios are extremely imbalanced by using the saddlepoint approximation to calibrate the distribution of score test statistics.

### Covariates

For each GWAS conducted for each phenotype and ancestry group, we included the following covariates:
- Age
- Sex
- Age * Sex
- Age<sup>2</sup>
- Age<sup>2</sup> * Sex
- The first 10 PCs
