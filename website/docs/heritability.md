---
id: heritability
title: Heritability Estimation
---

:::note
Important caveats with our heritability estimation approach can be found [here](https://pan.ukbb.broadinstitute.org/docs/qc#heritability).
:::note
:::note
Due to high computational cost associated with running genotype-based heritability estimation on the EUR ancestry group, we only provide genotype-based heritability estimates for a [pilot](ADDMANIFEST) set of traits in EUR, with results available for all non-EUR ancestry-trait pairs. Summary statistic methods were run for all ancestry-trait pairs.
:::note

We have used several methods to estimate heritability across up to 6 ancestry groups and 7,228 traits, totaling over 16,000 ancestry-trait pairs. All results can be found in the [phenotype manifest](ADDMANIFEST) and [heritability manifest](ADDMANIFEST). The following approaches were used:

- Univariate LD score regression ([LDSC](https://github.com/bulik/ldsc); `ldsc`), run using [LD score flat files](https://pan.ukbb.broadinstitute.org/docs/ld#ld-scores) using high-quality HapMap3 SNPs with MAF $\geq$ 0.01 with summary statistics exported from the [results table](https://pan.ukbb.broadinstitute.org/docs/hail-format#results-schema). See Bulik-Sullivan et al. 2015 Nat Gen for more information on the method.
- Stratified LD score regression ([S-LDSC](https://github.com/bulik/ldsc); `sldsc_25bin`), run using the same summary statistcs as `ldsc` with LD scores generated from SNPs in 5 MAF and 5 LD score bins. See Finucane et al. 2015 Nat Gen for more information on the method.
- Randomized Haseman-Elston ([RHEmc](https://github.com/sriramlab/RHE-mc); `rhemc_8bin`) using genotype data with 2 MAF and 4 LD score bins using default settings. This was predominantly used to analyze non-EUR ancestry groups but includes a small set of estimates for traits in EUR. See Pazokitoroudi et al. 2020 Nat Comm for more information on the method.
- Randomized Haseman-Elston (`rhemc_25bin`) using genotype data with 5 MAF and 5 LD score bins using default settings. Only run for non-EUR ancestry groups.
- Randomized Haseman-Elston (`rhemc_25bin_50rv`) using genotype data with 5 MAF and 5 LD score bins using 50 random variables to reduce run-to-run variability at a slightly higher computational cost. Only run for non-EUR ancestry groups.

Multi-component analyses were performed using either a 5x5 or 2x4 grid of MAF and LD score bins. MAF bins were either `(0, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.4), (0.4, 0.5)` or `(0.0, 0.05), (0.05, 0.5)`; LD score bins were computed as quantiles of the LD score distribution across all SNPs. LD scores used in creation of SNP bins were exported for *all SNPs with MAC > 20*; LD scores used for LD score regression were filtered after generation as described [here](https://pan.ukbb.broadinstitute.org/docs/ld#ld-scores).

All randomized Haseman-Elston runs were performed using the same fixed effects covariates used in the main GWAS; namely the first 20 genotype PCs, $age$, $sex$, $age*sex$, $age^2$, and $age^2*sex$. In the case of sex-specific phenotypes, only the first 20 genotype PCs, $age$, and $age^2$ were included.

Genotype and sample filtering was performed for genotype-based heritability estimation. We filtered to SNPs outside the MHC region that were defined with MAF $\geq 0.01$ for which we did not observe significant deviation from Hardy-Weinberg equilibrium. Importantly, we filtered to SNPs that passed the above criteria *for all ancestry groups*, ensuring the same set of SNPs for heritability estimation across all ancestry groups. We included only unrelated samples for genotype-based heritability estimation.

We note that our S-LDSC and LDSC-based estimates using summary statistics for individuals of EUR ancestry were highly concordant with [prior estimates](https://nealelab.github.io/UKBB_ldsc/) in UKB for overlapping phenotypes, however we observed very poor power for detection of $h^2 > 0$ for non-EUR ancestry groups. Genotype-based Haseman-Elston regression at scale (RHEmc) showed improved power for heritability dection in non-EUR ancestry groups and showed good concordance with S-LDSC in EUR, as discussed in our [post](https://pan.ukbb.broadinstitute.org/blog/2022/04/11/h2-qc-updated-sumstats).

We used these heritability estimates, along with other important summary statistics, to build a [comprehensive QC pipeline](https://pan.ukbb.broadinstitute.org/docs/qc#quality-control-of-summary-statistics).