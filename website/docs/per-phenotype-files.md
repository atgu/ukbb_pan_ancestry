---
id: per-phenotype-files
title: Per-phenotype files
---

## Overview

The data are released in 7,221 flat files, one for each phenotype, and a corresponding tabix index file for each. These files are available on Dropbox (for large-scale analysis, we recommend using the [Hail format](Hail-format) files on Google Cloud).

The files are named with respect to their `trait_type`, `phenocode`, and a combination of `pheno_sex`, `coding`, or `modifier`. To find a specific phenotype, we suggest looking in the [phenotype manifest (Google Sheets)](https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit?usp=sharing) (available for download on [Dropbox](https://www.dropbox.com/s/18p4lj3finj11oh/phenotype_manifest.tsv.bgz?dl=0)). Search for your phenotype(s) of interest and use the paths indicated to download the summary statistics. A description of fields in the manifest can be found [here](#phenotype-manifest-file).

The [per-phenotype files](#per-phenotype-files) are summary statistics files containing meta-analyzed and single-ancestry GWAS results. We especially highlight the `low_confidence` fields, which includes some (non-exhaustive) basic quality control filters (see below). These files each have 28,987,534 variants, but note that not all populations will have data for each variant.

Finally, the [variant manifest file](#variant-manifest-file) includes information on each variant in the dataset and has the same number of rows as each per-phenotype file. We highlight the `high_quality` column which represents variants that are PASS variants in gnomAD and have consistent frequencies with each population in gnomAD (AFR, AMR, EAS, and EUR frequencies are within 2-fold or chi-squared p-value of the difference > 1e-6).

## Phenotype manifest file

[Pan-UK Biobank phenotype manifest (Google Sheets)](https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit?usp=sharing) (download on [Dropbox](https://www.dropbox.com/s/18p4lj3finj11oh/phenotype_manifest.tsv.bgz?dl=0))

### Phenotype ID fields

The first 5 fields are guaranteed to be unique.

`trait_type`: One of the following: continuous, biomarkers, prescriptions, icd10, phecode, categorical  

`phenocode`: The code for the phenotype (for continuous, biomarkers, and categorical traits, this corresponds to the field ID as described by UKB, e.g. [21001 for BMI](http://biobank.ctsu.ox.ac.uk/showcase/field.cgi?id=21001))

`pheno_sex`: Indicating whether the phenotype was run for both sexes (`pheno_sex`="both_sexes") or in just females (`pheno_sex`="females") or males (`pheno_sex`="males"). In 0.1, this is only differentiated for phecodes.

`coding`: For categorical variables, this corresponds to the coding that was used (e.g. [coding 2](http://biobank.ctsu.ox.ac.uk/showcase/coding.cgi?id=100434) for [field 1747](http://biobank.ctsu.ox.ac.uk/showcase/field.cgi?id=1747)). For all other `trait_type`s, this field is blank.

`modifier`: Refers to any miscellaneous downstream modifications of the phenotype (e.g. `irnt` for inverse-rank normal transformation). If the phenotype is updated, this field can be used to denote the update (e.g. the particular wave of COVID-19 data used).

`description`: A shorter description of the phenotype (for continuous, biomarkers, and categorical variables, corresponds to the Description on the [showcase](http://biobank.ctsu.ox.ac.uk/showcase/field.cgi?id=30600)). For phecodes, this is the "description" column in the [phecodes definition file](https://github.com/atgu/ukbb_pan_ancestry/blob/master/data/PHECODE_v1.2b1_INFO_20200109.txt).

`description_more`: A longer description of the phenotype (for continuous and categorical variables, corresponds to the Notes page on the [showcase](http://biobank.ctsu.ox.ac.uk/showcase/field.cgi?id=30600)).

`coding_description`: For categorical variables, a description of the particular coding that was used (the Meaning column on the showcase page for that [coding](http://biobank.ctsu.ox.ac.uk/showcase/coding.cgi?id=100434)).

`category`: A categorization of the phenotype. For continuous, biomarkers, and categorical traits, this corresponds to the Category at the top of the [showcase page](http://biobank.ctsu.ox.ac.uk/showcase/field.cgi?id=30600). For ICD codes, this corresponds to the Chapter of the ICD code; for phecodes, this is the "group" column in the [phecodes definition file](https://github.com/atgu/ukbb_pan_ancestry/blob/master/data/PHECODE_v1.2b1_INFO_20200109.txt); for prescriptions, this corresponds to a semi-manual categorization of prescription drugs.

### Case and ancestry fields
:::note
If a trait is quantitative (`trait_type` is "continuous" or "biomarkers"), all samples are considered to be "cases". Thus, the number of cases is equivalent to the number of samples.
:::note

`n_cases_full_cohort_both_sexes`: Number of cases across all ancestry groups, females and males combined. 

`n_cases_full_cohort_females`: Number of female cases across all ancestry groups.

`n_cases_full_cohort_males`: Number of male cases across all ancestry groups.

`pops`: Comma-delimited list of ancestry codes for which this phenotypes was GWASed.

`num_pops`: Number of ancestry groups for which this phenotype was GWASed.

### Population-specific fields
:::note
The variable `pop` is a placeholder for a 3-letter ancestry code. For example, `n_cases_AFR` is the number of cases with AFR ancestry. 
:::note
:::note
If a trait is quantitative (`trait_type` is "continuous" or "biomarkers"), all samples are considered to be "cases". Thus, the number of cases is equivalent to the number of samples.
:::note

`n_cases_{pop}`: Number of cases with `pop` ancestry.

`n_controls_{pop}`: Number of controls with `pop` ancestry.

`saige_heritability_{pop}`: The heritability as estimated by SAIGE: note that this is likely not well-calibrated for binary traits, or traits with high heritabilities. A second estimate of heritability from LD score regression is coming soon.

`lambda_gc_{pop}`: The genomic control (lambda GC) calculated from the summary statistics for `pop` with low-confidence statistics removed and only considering high-quality variants.

### File information
:::note
For each field in this section there also exists a field with the suffix `_tabix`, which contains the equivalent information for the tabix file. For instance, `filename_tabix` contains the name of the tabix file.
:::note

`filename`: Name of summary statistics file.

`dropbox_link`: Dropbox link to download summary statistics file.

`wget`: wget command to download summary statistics file.

`size_in_bytes`: Size of summary statistics file in bytes.

`md5_hex`: MD5 hexadecimal hash.

## Per-phenotype files

The per-phenotype files are `tsv.bgz` files are (b)gzipped: they can either be unzipped (`zcat file.tsv.bgz > file.txt`), or read natively in R (`read_delim(gzfile('file.tsv.bgz'), delim='\t')`) and Python (`gzip.open('file.tsv.bgz')`).

Depending on whether a phenotype is quantitative (`trait_type` is "continuous" or "biomarkers") or binary (`trait_type` is "prescriptions", "icd10", "phecode" or "categorical"), the number of columns will change due to  case/control-stratified statistics for binary phenotypes.

### Variant fields

`chr`: Chromosome of the variant.

`pos`: Position of the variant in GRCh37 coordinates.

`ref`: Reference allele on the forward strand.

`alt`: Alternate allele (not necessarily minor allele). Used as effect allele for GWAS.

### Meta-analysis fields

`af_meta`: Alternate allele frequency from meta-analysis across populations for which this phenotype was GWASed. **NOTE: This field only appears in files for quantitative phenotypes.**

`af_cases_meta`: Alternate allele frequency in cases from meta-analysis across populations for which this phenotype was GWASed. **NOTE: This field only appears in files for binary phenotypes.**

`af_controls_meta`: Alternate allele frequency in controls from meta-analysis across populations for which this phenotype was GWASed. **NOTE: This field only appears in files for binary phenotypes.**

`beta_meta`: Estimated effect size of alternate allele from meta-analysis across populations for which this phenotype was GWASed.

`se_meta`: Estimated standard error of `beta_meta`.

`pval_meta`: p-value of `beta_meta` significance test.

`pval_heterogeneity`: p-value from heterogeneity test of meta-analysis.

### Population-specific fields
:::note
The variable `pop` used in this section is a placeholder for a 3-letter ancestry code. For example, `af_AFR` is the alternate allele frequency for AFR samples included in the GWAS of this phenotype.
:::note 

:::note
An ancestry-specific column is only included in the file if a GWAS was run for that ancestry. For example, a trait that was only GWASed in AMR and CSA samples will only have the fields `af_AMR`, `af_CSA`, `beta_AMR`, `beta_CSA`, etc.
:::note

`af_{pop}`: Alternate allele frequency for `pop` samples included in the GWAS of this phenotype. **NOTE: This field only appears in files for quantitative phenotypes.**

`af_cases_{pop}`: Alternate allele frequency for `pop` cases included in the GWAS of this phenotype. **NOTE: This field only appears in files for binary phenotypes.**

`af_controls_{pop}`: Alternate allele frequency for `pop` controls included in the GWAS of this phenotype. **NOTE: This field only appears in files for binary phenotypes.**


`beta_{pop}`: Estimated effect size of alternate allele from GWAS of `pop` samples.

`se_{pop}`: Estimated standard error of `beta_{pop}`.

`pval_{pop}`: p-value of `beta_{pop}` significance test.

`low_confidence_{pop}`: Boolean flag indicating low confidence for `pop` based on the following heuristics:
- Alternate allele count in cases <= 3
- Alternate allele count in controls <= 3
- Minor allele count (cases and controls combined) <= 20

## Variant manifest file

### Variant fields

As in [per-phenotype files](#Per-phenotype-files).

`chr`: Chromosome of the variant.

`pos`: Position of the variant in GRCh37 coordinates.

`ref`: Reference allele on the forward strand.

`alt`: Alternate allele (not necessarily minor allele). Used as effect allele for GWAS.

`rsid`: The RSID for the variant (from the BGEN file from UK Biobank).

`varid`: The variant ID for the variant (from the BGEN file from UK Biobank).

`pass_gnomad_genomes`: A boolean corresponding to the PASS status in gnomAD (`NA` if variant is not in gnomAD).

`n_passing_populations`: The number of populations (max 4: AFR, AMR, EAS, EUR) where the frequency in UKB is less than twice the frequency in gnomAD for the corresponding population (see below), and the p-value of a chi-squared test assessing the difference is > 1e-6.

`high_quality`: A boolean corresponding to a high-quality variant based on these filters (`pass_gnomad_genomes & n_passing_populations == 4`).

`nearest_genes`: The nearest genes for this variant based on Gencode v19.

`info`: The Info score for this variant (from the `ukb_mfi_chrN_v2.txt` files).

### Population-specific fields
:::note
The variable `pop` used in this section is a placeholder for a 3-letter ancestry code. For example, `af_AFR` is the alternate allele frequency for AFR individuals in the whole dataset.
:::note

`ac_{pop}`: The alternate allele count for this variant across all individuals in `pop`. Defined as `af_{pop} * an_{pop}`.

`af_{pop}`: The alternate allele frequency for this variant across all individuals in `pop`. The mean `dosage` divided by two.

`an_{pop}`: The alternate allele number for this variant. This is twice the number of `pop` individuals with a defined genotype at this site.

`gnomad_genomes_ac_{pop}`: The alternate allele count for this variant in the nearest gnomAD population: AFR, EAS, and AMR are matched as-is, while EUR is matched to the "North-West European" subset of gnomAD (not available for CSA or MID, as these populations are not in gnomAD v2 genomes).

`gnomad_genomes_af_{pop}`: The alternate allele frequency for this variant in the nearest gnomAD population.

`gnomad_genomes_an_{pop}`: The alternate allele number for this variant in the nearest gnomAD population. This is twice the number of individuals with a defined genotype at this site.

