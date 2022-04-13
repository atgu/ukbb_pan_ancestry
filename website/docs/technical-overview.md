---
id: technical-overview
title: Overview
---


## Pan-ancestry GWAS of UK Biobank

Here, we present a multi-ancestry analysis of 7,221 phenotypes using a generalized mixed model association testing framework, spanning 16,119 genome-wide association studies. We provide standard meta-analysis across all populations and with a leave-one-population-out approach for each trait. We develop a stringent quality control pipeline, identifying variants that are discrepant with gnomAD frequencies, and make recommendations for filtering these and other GWAS results.

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

We release the summary statistics in two formats:

- For one or a few phenotypes, we recommend using the phenotype-specific flat files: see further description [here](per-phenotype-files).

- For analysis the full dataset (all phenotypes, all populations), the summary statistics are available in Hail formats: see further description [here](hail-format).

### Approach

Analysis was done using [SAIGE](https://github.com/weizhouUMICH/SAIGE/wiki/Genetic-association-tests-using-SAIGE) implemented in [Hail Batch](https://hail.is/docs/batch/index.html) to parallelize across populations, phenotypes, and regions of the genome. More details can be found below:

- Details about the QC process can be found [here](qc) including [determination of ancestry groups](qc#ancestry-definitions).
- Description of GWAS pipeline and implementation can be found on our [Github](https://github.com/atgu/ukbb_pan_ancestry/wiki/Batch-pipeline).

The sample size for each population and the number of phenotypes run is as follows:

| Population | Num. Individuals | Num. Phenotypes |
|-----|-----------|----------|
| AFR |      6636 |     2493 |
| AMR |       980 |     1105 |
| CSA |      8876 |     2771 |
| EAS |      2709 |     1612 |
| EUR |    420531 |     7200 |
| MID |      1599 |     1372 |

Each phenotype may have fewer samples run, depending on data missingness, which can be found in the [phenotype manifest](https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=903887429), or `n_cases` and `n_controls` in the Hail MatrixTable.
