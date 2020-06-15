---
id: technical-overview
title: Technical overview
---


## Pan-ancestry GWAS of UK Biobank

The UK Biobank is a collection of a half million individuals with paired genetic and phenotype information that has been enormously valuable in studies of genetic etiology for common diseases and traits. However, most genome-wide analyses of this dataset use only the European ancestry individuals in an effort to control for population structure. Analyzing a more inclusive and diverse dataset increases power and improves the potential for discovery. Further, the analysis of all ancestries offers the opportunity to evaluate similarity in genetic architecture as well as enhance the generalizability of genetic discoveries.

Here, we present a multi-ancestry analysis of 7,256 phenotypes using a generalized mixed model association testing framework, spanning 16,523 genome-wide association studies. We provide standard meta-analysis across all populations and with a leave-one-population-out approach for each trait. We develop a stringent quality control pipeline, identifying variants that are discrepant with gnomAD frequencies, and make recommendations for filtering these and other GWAS results. Here, we release these summary statistics freely to the community ahead of publication.

With knowledge that these analyses include some culturally sensitive phenotypes, we developed [FAQs]() as a resource to understand our study, its risks versus benefits, and limitations for interpreting results in collaboration with local affinity groups and a bioethics expert. We intend for this resource to promote more inclusive research practices among researchers, to accelerate novel scientific discoveries that are more generalizable, and improve the health of all people equitably.

### Multi-ancestry analysis

Participants have been divided into ancestry groups to account for population stratification in GWAS analyses. Throughout these docs, these ancestry groupings are referred to by 3-letter ancestry codes derived from or closely related to those used in the 1000 Genomes Project and Human Genome Diversity Panel, as follows:

- `EUR` = European ancestry
- `CSA` = Central/South Asian ancestry
- `AFR` = African ancestry
- `EAS` = East Asian ancestry
- `MID` = Middle Eastern ancestry
- `AMR` = Admixed American ancestry
 
These codes refer only to ancestry groupings used in GWAS, not necessarily other demographic or self-reported data.

### Release data

We are planning to release the summary statistics in two formats:

- For one or a few phenotypes, we recommend using the phenotype-specific flat files: see further description [here](https://github.com/atgu/ukbb_pan_ancestry/wiki/Per-phenotype-files).

- For analysis the full dataset (all phenotypes, all populations), the summary statistics are available in Hail formats: see further description [here](https://github.com/atgu/ukbb_pan_ancestry/wiki/Hail-format).

### Approach

Analysis was done using [SAIGE](https://github.com/weizhouUMICH/SAIGE/wiki/Genetic-association-tests-using-SAIGE) implemented in [Hail batch](https://hail.is/docs/batch/index.html) to parallelize across populations, phenotypes, and regions of the genome. More details can be found below:

- Details about the QC process can be found [here](https://github.com/atgu/ukbb_pan_ancestry/wiki/QC) including [determination of ancestry groups](https://github.com/atgu/ukbb_pan_ancestry/wiki/QC#ancestry-definitions).
- Description of GWAS pipeline and implementation is [here](https://github.com/atgu/ukbb_pan_ancestry/wiki/Batch-pipeline).

The sample size for each population and the number of phenotypes run is as follows:

| Population | Num. Individuals | Num. Phenotypes |
|-----|-----------|----------|
| AFR |      6636 |     2484 |
| AMR |       980 |     1105 |
| CSA |      8876 |     2765 |
| EAS |      2709 |     1612 |
| EUR |    420531 |     7185 |
| MID |      1599 |     1372 |

Each phenotype may have fewer samples run, depending on data missingness, which can be found in the [phenotype manifest](), or `n_cases` and `n_controls` in the Hail MatrixTable.
