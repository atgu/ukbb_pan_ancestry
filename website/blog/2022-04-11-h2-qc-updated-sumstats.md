---
id: h2-qc-updated-sumstats
title: Quality control, heritability analyses, and updates to summary statistics
author: Rahul Gupta, on behalf of the Pan UKBB Team
tags: [sumstats, heritability, independent phenotypes, meta analysis]
---

We are excited to report significant updates to our summary statistics and data release:

1. We performed heritability analyses across > 16,000 ancestry-trait pairs using several approaches.
2. We developed a detailed summary statistics QC approach to prioritize the highest-quality phenotypes best suited for downstream analyses.
3. We identified a maximally independent set of phenotypes that passed our QC filters.
4. We recomputed summary statistics for traits that showed extremely significant p-values with standard errors of 0, now with non-zero standard errors and $\log p$-values to avoid numerical underflow.
5. We updated cross-ancestry meta-analyses to incorporate updated summary statistics and also computed new meta-analyses using only QC-pass ancestry-trait pairs.

<!--truncate-->

## How to use QC-ed sumstats

If using the UKBB Pan Ancestry [codebase](https://github.com/atgu/ukbb_pan_ancestry), the final summary statistics MatrixTable filtered to the phenotype-ancestry pairs passing all QC can be obtained via:

```
mt = load_final_sumstats_mt(filter_pheno_h2_qc=True)
```

If using the manifest flat file directly (available [here](https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=903887429)), filter on `phenotype_qc_{ancestry} == PASS` for the ancestry of interest.

## Heritability analyses

We produced narrow-sense heritability estimates using univariate linkage-disequilibrium score regression (LDSC) leveraging LD scores computed as [previously described](https://pan.ukbb.broadinstitute.org/blog/2020/10/29/ld-release), as well as using stratified LDSC with 25 minor allele frequency (MAF) and LD score bins. Our results for summary statistics for individuals of EUR ancestry were highly concordant with [prior estimates](https://nealelab.github.io/UKBB_ldsc/) in UKB for overlapping phenotypes (Figure 1), however we observed very poor power for detection of $h^2 > 0$ for non-EUR ancestry groups.

<center><img src="/img/EUR_SLDSC_LDSC_compare_all_phenotypes.png"  width="608" height="600"/></center>

**Figure 1**: Comparison of liability-scale $h^2$ estimates for EUR ancestry-trait pairs using S-LDSC (25 bins) in UKB versus prior estimates using S-LDSC with the baselineLDv1.1 model in UKB among EUR samples. Black line is $y=x$.

To boost power for non-EUR ancestry groups we leveraged Haseman-Elston regression at scale on genotypes implemented in RHEmc ([Pazokitoroudi et al. 2020 Nat Comm](https://www.nature.com/articles/s41467-020-17576-9)). We piloted this approach for several [pilot phenotypes](https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=1101141753) expected to behave well, finding good correlations between estimates in EUR using RHEmc and S-LDSC (Figure 2).

<center><img src="/img/EUR_rhemc_nonliab_v_sldsc_allpananc.png"  width="640" height="523"/></center>

**Figure 2**: Observed-scale $h^2$ estimates in select phenotypes among EUR individuals obtained using RHEmc (8 bins, to limit computational cost) vs. using S-LDSC (25 bins), colored by trait type category. Black line is $y=x$.

We then extended this approach to all non-EUR ancestry-trait pairs and report these results using both 8 and 25 MAF-LD score bins. We observe an increase in power for detecting non-zero $h^2$ using this method for non-EUR individuals relative to S-LDSC (Figure 3) as expected from genotype-based approaches. Due to computational cost and the similarity of $h^2$ estimates observed for EUR individuals between S-LDSC and RHEmc, our final estimates use S-LDSC for EUR and RHEmc for non-EUR ancestry groups.

<center><img src="/img/h2z_densities_nonEUR_sldsc_rhemc.png"  width="590" height="413"/></center>

**Figure 3**: Distribution of $h^2$ z-scores for non-EUR ancestry-trait pairs in UKB when computed using summary statistics (S-LDSC, pink) or genotype data (RHEmc, blue). Dotted line corresponds to $Z = 4$.

We have made our final heritability estimates available for download in the [phenotype manifest](https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=903887429) and data from all methods available in the [heritability manifest](https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=1797288938). All results are also available in [MatrixTable format](https://pan.ukbb.broadinstitute.org/downloads). Further description of the heritability estimation approach can be found [here](https://pan.ukbb.broadinstitute.org/docs/heritability).

## QC approach

We observed significantly negative RHEmc heritability estimates for traits that were likely to have a high propensity of confounding (e.g., UK geographic location). We subsequently developed a systematic QC approach using the following *sequential* filters to prioritize high quality phenotypes for further analysis:

1. `GWAS_run`: if the GWAS was performed for the ancestry-trait pair
2. `ancestry_reasonable_n`: if the ancestry has a reasonable sample size; has the effect of removing AMR
3. `defined_h2`: if the heritability estimate is non-missing
4. `significant_z`: if the ancestry-trait pair shows $h^2$ z-score $> 0$
5. `in_bounds_h2`: if, for all ancestries for a given trait, observed-scale heritability estimates $\in (0,1)$
6. `normal_lambda`: if, for all ancestries for a given trait, $\lambda_{GC} > 0.9$
7. `normal_ratio`: if, for the top three best powered ancestry groups (EUR, CSA, AFR), the S-LDSC ratio, given by $\frac{intercept-1}{mean \chi^2 -1}$, $< 0.3$ or the ratio z-score $< 4$
8. `EUR_plus_1`: if the trait passes all above filters in EUR and at least 1 other ancestry group

A more detailed description of the approach is forthcoming. All QC results for ancestry-trait pairs can be found in the [phenotype manifest](https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=903887429) and [heritability manifest](https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=1797288938).

## Maximally independent set

Using the set of phenotypes with ancestries passing all of the above QC filters, we next wanted to construct a set of independent phenotypes to avoid double-counting in downstream analyses. To this end we constructed a maximally independent set of phenotypes using a pairwise phenotypic correlation matrix (released [here](https://pan.ukbb.broadinstitute.org/downloads) via [`make_pairwise_ht`](https://github.com/Nealelab/ukb_common/blob/f9b4c037b57e932a52dcfb8c35f1e077c6939610/src/ukbb_common/utils/phenotype_loading.py#L338)). Of all phenotype pairs, we retained any with a pairwise correlation $r < 0.1$. For phenotype pairs with  $r > 0.1$ , we used [`hl.maximal_independent_set`](https://hail.is/docs/0.2/methods/misc.html#hail.methods.maximal_independent_set) to identify indendent phenotypes for retention, imposing a tiebreaker of higher case count (or higher sample size for continuous phenotypes). This resulted in the selection of 195 independent phenotypes, which were used in subsequent analyses analyzing trends across phenotypes. We have released this set in the [phenotype manifest](https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=903887429).

## Updated summary statistics and meta analyses

Using the updated summary statistics (now with $\log p$-values), we reproduced cross-ancestry meta-analyses for all phenotypes (using only ancestry groups for which GWAS was completed). These are represented in the corresponding phenotype-specific summary statistics as columns with the `_meta` suffix.

We also produced a new set of "high quality" meta-analyses which only include ancestries passing QC (see **QC approach**). These meta-analyses are only available for phenotypes for which QC-passing ancestry groups exist. These new "high quality" meta-analysis results are present in columns with the `_meta_hq` suffix in the corresponding summary statistics for phenotypes for which they were computed.

We have updated our [codebase](https://github.com/atgu/ukbb_pan_ancestry) to export these updated manifests and summary statistics, which have now been uploaded to [AWS](https://pan.ukbb.broadinstitute.org/downloads). See the [summary statistics key](https://pan.ukbb.broadinstitute.org/docs/per-phenotype-files) for information on the summary statistics columns containing these meta-analyses as well as other new fields.

Please note that the previous iteration of release files have been archived at `s3://pan-ukb-us-east-1/archive_20200615/`.