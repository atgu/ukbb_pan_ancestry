---
id: per-phenotype-files
title: Per-phenotype files
---

## Overview

The data are released in 7,228 flat files, one for each phenotype, and a corresponding tabix index file for each. These files are available on Amazon AWS (for large-scale analysis, we recommend using the [Hail format](Hail-format) files on Google Cloud).

The files are named with respect to their `trait_type`, `phenocode`, and a combination of `pheno_sex`, `coding`, or `modifier`. To find a specific phenotype, we suggest looking in the [phenotype manifest (Google Sheets)](ADDMANIFEST) (available for download on [Amazon](https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/phenotype_manifest.tsv.bgz)). Search for your phenotype(s) of interest and use the paths indicated to download the summary statistics. A description of fields in the manifest can be found [here](#phenotype-manifest-file).

The [per-phenotype files](#per-phenotype-files) are summary statistics files containing meta-analyzed and single-ancestry GWAS results. We especially highlight the `low_confidence` fields, which includes some (non-exhaustive) basic quality control filters (see below). These files each have 28,987,534 variants, but note that not all populations will have data for each variant.

Finally, the [variant manifest file](#variant-manifest-file) includes information on each variant in the dataset and has the same number of rows as each per-phenotype file. We highlight the `high_quality` column which represents variants that are PASS variants in gnomAD and have consistent frequencies with each population in gnomAD (AFR, AMR, EAS, and EUR frequencies are within 2-fold or chi-squared p-value of the difference > 1e-6).

## Phenotype manifest file

[Pan-UK Biobank phenotype manifest (Google Sheets)](ADDMANIFEST) (download on [Amazon](https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/phenotype_manifest.tsv.bgz))

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

`in_max_independent_set`: If the phenotype is in our [maximally indepdent set](https://pan.ukbb.broadinstitute.org/blog/2022/04/11/h2-qc-updated-sumstats). This set of relatively uncorrelated phenotypes was constructed using a pairwise phenotypic correlation matrix of phenotypes with ancestries passing all QC filters (released [here](https://pan.ukbb.broadinstitute.org/downloads) via [`make_pairwise_ht`](https://github.com/Nealelab/ukb_common/blob/f9b4c037b57e932a52dcfb8c35f1e077c6939610/src/ukbb_common/utils/phenotype_loading.py#L338)). Of all phenotype pairs, we retained any with a pairwise correlation $r < 0.1$. For pairs with  $r > 0.1$ , we used [`hl.maximal_independent_set`](https://hail.is/docs/0.2/methods/misc.html#hail.methods.maximal_independent_set) to identify indendent phenotypes for retention, imposing a tiebreaker of higher case count (or higher sample size for continuous phenotypes), producing 195 independent phenotypes.

### Case and ancestry fields
:::note
If a trait is quantitative (`trait_type` is "continuous" or "biomarkers"), all samples are considered to be "cases". Thus, the number of cases is equivalent to the number of samples.
:::note

`n_cases_full_cohort_both_sexes`: Number of cases (or individuals phenotyped for quantitative traits) across all ancestry groups, females and males combined. Should be similar to the sum of per-ancestry n_cases for relevant ancestries, but may include ancestry outliers and samples that failed QC.

`n_cases_full_cohort_females`: Number of female cases (or individuals phenotyped for quantitative traits) across all ancestry groups. May include ancestry outliers and samples that failed QC.

`n_cases_full_cohort_males`: Number of male cases (or individuals phenotyped for quantitative traits) across all ancestry groups. May include ancestry outliers and samples that failed QC.

`n_cases_hq_cohort_both_sexes `: Number of cases (or individuals phenotyped for quantitative traits) across ancestry groups passing stringent phenotype QC (see `pops_pass_qc`), females and males combined. Should be similar to the sum of per-ancestry n_cases for relevant ancestries, but may include ancestry outliers and samples that failed QC.

`n_cases_hq_cohort_females`: Number of female cases (or individuals phenotyped for quantitative traits) across ancestry groups passing stringent phenotype QC (see `pops_pass_qc`). May include ancestry outliers and samples that failed QC.

`n_cases_hq_cohort_males`: Number of male cases (or individuals phenotyped for quantitative traits) across ancestry groups passing stringent phenotype QC (see `pops_pass_qc`). May include ancestry outliers and samples that failed QC.

`pops`: Comma-delimited list of ancestry codes for which this phenotypes was GWASed.

`num_pops`: Number of ancestry groups for which this phenotype was GWASed.

`pops_pass_qc`: Comma-delimited list of ancestry codes for which this phenotype passes QC (see [quality control](https://pan.ukbb.broadinstitute.org/docs/qc#quality-control-of-summary-statistics), [heritability manifest](ADDMANIFEST), and `phenotype_qc_{pop}` field).

`num_pops_pass_qc`: Number of ancestry groups for which this phenotype passes QC.

### Population-specific fields
:::note
The variable `pop` is a placeholder for a 3-letter ancestry code. For example, `n_cases_AFR` is the number of cases with AFR ancestry. 
:::note
:::note
If a trait is quantitative (`trait_type` is "continuous" or "biomarkers"), all samples are considered to be "cases". Thus, the number of cases is equivalent to the number of samples.
:::note
:::note
The final heritability estimates for non-EUR traits use a randomized Haseman Elston estimator on genotype data with 5 MAF and 5 LD bins and 50 random variables (`rhemc_25bin_50rv`), while for computational tractability for EUR we use S-LDSC with 5 MAF and 5 LD bins. In the phenotype manifest we provide estimates with the `sldsc_25bin_` prefix for EUR and the `rhemc_25bin_50rv_` prefix for other ancestry groups. See the [heritability manifest](ADDMANIFEST) for all heritability estimates computed, including `sldsc_25bin` for all ancestry groups.
:::note

`n_cases_{pop}`: Number of cases (or individuals phenotyped for quantitative traits) with `pop` ancestry in the GWAS analysis. Excludes ancestry outliers and samples that failed QC.

`n_controls_{pop}`: Number of controls with `pop` ancestry in the GWAS analysis. Excludes ancestry outliers and samples that failed QC.

`{rhemc_25bin_50rv/sldsc_25bin}_h2_observed_{pop}`: Observed scale heritability estimates using 25 MAF/LD bins with RHEmc (non-EUR) or SLDSC (EUR).

`{rhemc_25bin_50rv/sldsc_25bin}_h2_observed_se_{pop}`: Observed scale heritability standard error using 25 MAF/LD bins with RHEmc (non-EUR) or SLDSC (EUR).

`{rhemc_25bin_50rv/sldsc_25bin}_h2_liability_{pop}`: Libaility scale heritability estimates using 25 MAF/LD bins with RHEmc (non-EUR) or SLDSC (EUR), transformed using per-ancestry in-sample prevalence.

`{rhemc_25bin_50rv/sldsc_25bin}_h2_liability_se_{pop}`: Liability scale heritability standard error using 25 MAF/LD bins with RHEmc (non-EUR) or SLDSC (EUR), transformed using per-ancestry in-sample prevalence.

`{rhemc_25bin_50rv/sldsc_25bin}_h2_z_{pop}`: Heritability Z-scores (for test of h2 > 0) using per-ancestry-trait pair h2 estimates and standard errors.

`lambda_gc_{pop}`: The genomic control (lambda GC) calculated from the summary statistics for `pop` with low-confidence statistics removed and only considering high-quality variants.

`phenotype_qc_{pop}`: Phenotype QC outcome for each ancestry-trait pair. Filters are described in the [heritability manifest](ADDMANIFEST) in more detail. Filters are applied sequentially; this field specifies either PASS or the reason for failure.

### File information
:::note
For each field in this section there also exists a field with the suffix `_tabix`, which contains the equivalent information for the tabix file. For instance, `filename_tabix` contains the name of the tabix file.
:::note

`filename`: Name of summary statistics file.

`aws_link`: Link to download summary statistics file from Amazon AWS.

## Heritability manifest file

[Pan-UK Biobank heritability manifest (Google Sheets)](ADDMANIFEST) (download on [Amazon](https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/h2_manifest.tsv.bgz))

### Phenotype ID fields

:::note
The variable `pop` is a placeholder for a 3-letter ancestry code. For example, `n_cases_AFR` is the number of cases with AFR ancestry. 
:::note
:::note
Unlike the phenotype manifest, the heritability manifest provides a row for each ancestry group. Thus the *first 6 fields* are guarenteed to be unique.
:::note

`trait_type`: One of the following: continuous, biomarkers, prescriptions, icd10, phecode, categorical  

`phenocode`: The code for the phenotype (for continuous, biomarkers, and categorical traits, this corresponds to the field ID as described by UKB, e.g. [21001 for BMI](http://biobank.ctsu.ox.ac.uk/showcase/field.cgi?id=21001))

`pheno_sex`: Indicating whether the phenotype was run for both sexes (`pheno_sex`="both_sexes") or in just females (`pheno_sex`="females") or males (`pheno_sex`="males"). In 0.1, this is only differentiated for phecodes.

`coding`: For categorical variables, this corresponds to the coding that was used (e.g. [coding 2](http://biobank.ctsu.ox.ac.uk/showcase/coding.cgi?id=100434) for [field 1747](http://biobank.ctsu.ox.ac.uk/showcase/field.cgi?id=1747)). For all other `trait_type`s, this field is blank.

`modifier`: Refers to any miscellaneous downstream modifications of the phenotype (e.g. `irnt` for inverse-rank normal transformation). If the phenotype is updated, this field can be used to denote the update (e.g. the particular wave of COVID-19 data used).

`heritability.pop`: Ancestry group.

### Heritability methods

:::note
Each of the below **Heritability methods** subheadings contain (a subset of) the fields described in the **Heritability estimates** section (except for SAIGE estimates, which are provided as-is).
:::note
:::note
`rhemc_8bin` was run only for a subset of traits in the EUR ancestry group and for all traits across non-EUR ancestry groups. `rhemc_25bin` and `rhemc_25bin_50rv` were run only for non-EUR ancestry groups due to computational complexity. `sldsc_25bin` and `ldsc` were run for all ancestry-trait pairs.
:::note
:::note
More information can be found [here](https://pan.ukbb.broadinstitute.org/docs/heritability) on the heritability estimation approach, and important caveats can be found [here](https://pan.ukbb.broadinstitute.org/docs/qc#heritability).
:::note

`heritability.estimates.ldsc.*`: Univariate LD score regression

`heritability.estimates.sldsc_25bin.*`: Stratified LD score regression, 5 LD score bins x 5 MAF bins

`heritability.estimates.rhemc_25bin.*`: RHEmc (HE regression), 5 LD score bins x 5 MAF bins

`heritability.estimates.rhemc_8bin.*`: RHEmc (HE regression), 4 LD score bins x 2 MAF bins

`heritability.estimates.rhemc_25bin_50rv.*`:  RHEmc (HE regression), 5 LD score bins x 5 MAF bins; 50 random variables for improved power

`heritability.estimates.saige`: The heritability as estimated by SAIGE: note that this is likely not well-calibrated for binary traits, or traits with high heritabilities.

`heritability.estimates.final.*`: Final estimates; 25 bin SLDSC for EUR and 25 bin, 50 RV RHEmc for non-EUR (these are also present in the full manifest)

### Heritability estimates

`heritability.estimates.*.h2_observed`: Observed scale heritability point estimates.

`heritability.estimates.*.h2_liability`: Liability scale results, using per-ancestry in-sample prevalence.

`heritability.estimates.*.h2_z`: Heritability Z-scores for test of $h^2 > 0$.

`heritability.estimates.*.intercept`: LDSC intercept. Only present for `ldsc` and `sldsc_25bin`

`heritability.estimates.*.ratio`: LDSC ratio, given by $\frac{intercept-1}{mean \chi^2 -1}$ as a measure of the proportion of test statistic inflation not attributed to polygenicity. Only present for `ldsc` and `sldsc_25bin`.

### QC flags

:::note
The below QC flags were applied sequentially. See [quality control](https://pan.ukbb.broadinstitute.org/docs/qc#quality-control-of-summary-statistics) for more information.
:::note

`heritability.N_ancestry_QC_pass`: Number of ancestries passing all QC per trait.

`heritability.qcflags.GWAS_run`: if the GWAS was performed for the ancestry-trait pair

`heritability.qcflags.ancestry_reasonable_n`: if the ancestry has a reasonable sample size; has the effect of removing AMR

`heritability.qcflags.defined_h2`: if the heritability estimate is non-missing

`heritability.qcflags.significant_z`: if the ancestry-trait pair shows $h^2$ z-score $> 0$

`heritability.qcflags.in_bounds_h2`: if, for all ancestries for a given trait, observed-scale heritability estimates $\in (0,1)$

`heritability.qcflags.normal_lambda`: if, for all ancestries for a given trait, $\lambda_{GC} > 0.9$

`heritability.qcflags.normal_ratio`: if, for the top three best powered ancestry groups (EUR, CSA, AFR), the S-LDSC ratio, given by $\frac{intercept-1}{mean \chi^2 -1}$, $< 0.3$ or the ratio z-score $< 4$

`heritability.qcflags.EUR_plus_1`: if the trait passes all above filters in EUR and at least 1 other ancestry group

`heritability.qcflags.pass_all`: if all QC flags pass for each ancestry-trait pair

## Per-phenotype files

The per-phenotype files are `tsv.bgz` files are (b)gzipped: they can either be unzipped (`zcat file.tsv.bgz > file.txt`), or read natively in R (`read_delim(gzfile('file.tsv.bgz'), delim='\t')`) and Python (`gzip.open('file.tsv.bgz')`).

Depending on whether a phenotype is quantitative (`trait_type` is "continuous" or "biomarkers") or binary (`trait_type` is "prescriptions", "icd10", "phecode" or "categorical"), the number of columns will change due to  case/control-stratified statistics for binary phenotypes.

### Variant fields

`chr`: Chromosome of the variant.

`pos`: Position of the variant in GRCh37 coordinates.

`ref`: Reference allele on the forward strand.

`alt`: Alternate allele (not necessarily minor allele). Used as effect allele for GWAS.

### Meta-analysis fields

:::note
All meta-analyses were only performed on variants that were not flagged as `low_confidence` for a given population. As described below in **Population-specific fields**, per-variant `low_confidence` status is specific to each ancestry-trait pair.
:::note

`af_meta`: Alternate allele frequency from meta-analysis across populations for which this phenotype was GWASed. **NOTE: This field only appears in files for quantitative phenotypes.**

`af_cases_meta`: Alternate allele frequency in cases from meta-analysis across populations for which this phenotype was GWASed. **NOTE: This field only appears in files for binary phenotypes.**

`af_controls_meta`: Alternate allele frequency in controls from meta-analysis across populations for which this phenotype was GWASed. **NOTE: This field only appears in files for binary phenotypes.**

`beta_meta`: Estimated effect size of alternate allele from meta-analysis across populations for which this phenotype was GWASed.

`se_meta`: Estimated standard error of `beta_meta`.

`pval_meta`: p-value of `beta_meta` significance test.

`pval_heterogeneity`: p-value from heterogeneity test of meta-analysis.

### High quality meta-analysis fields

:::note
These fields are only present in flat files for traits that have any QC-pass ancestries. As a requirement for passing QC a trait must pass in EUR and at least 1 other ancestry (see [quality control](https://pan.ukbb.broadinstitute.org/docs/qc#quality-control-of-summary-statistics)).
:::note
:::note
As above, meta-analyses were only performed on variants that were not flagged as `low_confidence` for a given population. As described below in **Population-specific fields**, per-variant `low_confidence` status is specific to each ancestry-trait pair.
:::note

`af_meta_hq`: Alternate allele frequency from meta-analysis across populations for which this phenotype passes all QC filters. **NOTE: This field only appears in files for quantitative phenotypes.**

`af_cases_meta_hq`: Alternate allele frequency in cases from meta-analysis across populations for which this phenotype passes all QC filters. **NOTE: This field only appears in files for binary phenotypes.**

`af_controls_meta_hq`: Alternate allele frequency in controls from meta-analysis across populations for which this phenotype passes all QC filters. **NOTE: This field only appears in files for binary phenotypes.**

`beta_meta_hq`: Estimated effect size of alternate allele from meta-analysis across populations for which this phenotype passes all QC filters.

`se_meta_hq`: Estimated standard error of `beta_meta_hq`.

`pval_meta_hq`: p-value of `beta_meta_hq` significance test.

`pval_heterogeneity_hq`: p-value from heterogeneity test of meta-analysis.

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

Variant manifest (download from [Amazon AWS](https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/full_variant_qc_metrics.txt.bgz))

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

