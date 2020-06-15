---
id: per-phenotype-files
title: Per-phenotype files
---

[Per-phenotype files](#Per-phenotype-file-fields) are summary statistics files containing meta-analyzed and single-ancestry GWAS results. The files are in: `gs://ukb-diverse-pops/sumstats_flat_files/`

Phenotype information can be found in the phenotype manifest: `gs://ukb-diverse-pops/combined_results/phenotype_manifest.tsv.bgz`

## Phenotype manifest fields

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

`n_cases_full_cohort_both_sexes`: Number of cases across all ancestry groups, females and males combined.

`n_cases_full_cohort_females`: Number of female cases across all ancestry groups.

`n_cases_full_cohort_males`: Number of male cases across all ancestry groups.

`pops`: Comma-delimited list of ancestry codes for which this phenotypes was GWASed.

`num_pops`: Number of ancestry groups for which this phenotype was GWASed.

### Population-specific fields
> Note: The variable `pop` is a placeholder for a 3-letter ancestry code. For example, `n_cases_AFR` is the number of cases with AFR ancestry. 

`n_cases_{pop}`: Number of cases with `pop` ancestry.

`n_controls_{pop}`: Number of controls with `pop` ancestry.

`saige_heritability_{pop}`: The heritability as estimated by SAIGE: note that this is likely not well-calibrated for binary traits, or traits with high heritabilities. A second estimate of heritability from LD score regression is coming soon.

`lambda_gc_{pop}`: The genomic control (lambda GC) calculated from the summary statistics for `pop` with low-confidence statistics removed and only considering high-quality variants.

### File location

`filename`: Name of summary statistics file. All files are block compressed and tab delimited.

`dropbox_link`: Location of file as Dropbox link.

# Per-phenotype file fields

### Variant fields

`chr`: Chromosome of the variant.

`pos`: Position of the variant in GRCh37 coordinates.

`ref`: Reference allele on the forward strand.

`alt`: Alternate allele (not necessarily minor allele). Used as effect allele for GWAS.

### Meta-analysis fields

`af_meta`: Alternate allele frequency from meta-analysis across populations for which this phenotype was GWASed.

`beta_meta`: Estimated effect size of alternate allele from meta-analysis across populations for which this phenotype was GWASed.

`se_meta`: Estimated standard error of `beta_meta`.

`pval_meta`: p-value of `beta_meta` significance test.

`pval_heterogeneity`: p-value from heterogeneity test of meta-analysis.

### Population-specific fields
> Note: The variable `pop` used in this section is a placeholder for a 3-letter ancestry code. For example, `af_AFR` is the alternate allele frequency for AFR samples included in the GWAS of this phenotype. 

> Note: An ancestry-specific column is only included in the file if a GWAS was run for that ancestry. For example, a trait that was only GWASed in AMR and CSA samples will only have the fields `af_AMR`, `af_CSA`, `beta_AMR`, `beta_CSA`, etc.

`af_{pop}`: Alternate allele frequency for `pop` samples included in the GWAS of this phenotype.

`beta_{pop}`: Estimated effect size of alternate allele from GWAS of `pop` samples.

`se_{pop}`: Estimated standard error of `beta_{pop}`.

`pval_{pop}`: p-value of `beta_{pop}` significance test.

`low_confidence_{pop}`: Boolean flag indicating low confidence for `pop` based on the following heuristics:
- Alternate allele count in cases <= 3
- Alternate allele count in controls <= 3
- Minor allele count (cases and controls combined) <= 20


