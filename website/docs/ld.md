---
id: ld
title: LD scores and matrices
---

## Overview

We computed in-sample dosage-based LD matrices and scores for each of six ancestry group in UKBB. [LD matrices](https://pan.ukbb.broadinstitute.org/docs/ld#ld-matrices) are available in Hail's [BlockMatrix](https://hail.is/docs/0.2/linalg/hail.linalg.BlockMatrix.html) format on Amazon AWS (see details [here](http://localhost:3000/docs/hail-format)). [LD scores](https://pan.ukbb.broadinstitute.org/docs/ld#ld-scores) are available in [LDSC](https://github.com/bulik/ldsc)-compatible flat files (`.l2.ldscore.gz` and `.M_5_50`) [here](https://pan-ukb-us-east-1.s3.amazonaws.com/ld_release/UKBB.ALL.ldscore.tar.gz). For large-scale analysis, you can also find a full LD score [HailTable](https://hail.is/docs/0.2/hail.Table.html) (not restricted to the HapMap3 variants) on Amazon AWS (see details [here](http://localhost:3000/docs/hail-format))

For LD computation, please find technical details below. All the code is also publicly available [here](https://github.com/atgu/ukbb_pan_ancestry/blob/master/compute_ld_matrix.py). Detailed instruction for how to run LD score regression is available on [LDSC's website](https://github.com/bulik/ldsc/wiki).

## LD matrices

* The dosage-based genotype matrix $X$ was column-wise mean-centered and normalized.
* We applied the same variant QC filter used for the Pan-UKB GWAS (INFO > 0.8, MAC > 20 in each population; see details [here](https://pan.ukbb.broadinstitute.org/docs/qc#variant-qc))
* For covariate correction, the residuals from the regression of $genotype \sim covariates$ were obtained via $X_{adj} = M_cX$ where $M_c = I - C(C^TC)^{-1}C^T$, the residual-maker matrix, and $C$ is the matrix of covariates.
* We used the same covariates used for the Pan-UKB GWAS, namely $age$, $sex$, $age*sex$, $age^2$, $age^2*sex$, and the first 10 PCs of the genotype matrix (see details [here](https://pan.ukbb.broadinstitute.org/docs/qc#gwas-model)).
* We then computed LD matrix $R$ via $R = \frac{X_{adj}^TX_{adj}}{n}$ with a radius of <b><u>10 Mb</u></b>. Each element $\hat{r}_{jk}$ of $R$ represents the Pearson correlation coefficient of genotypes between variant $j$ and $k$.
* For X-chromosome, we computed a LD matrix jointly using both males and females where male genotypes are coded 0/1 and female genotypes are coded 0/1/2.

## LD scores

* To account for an upward bias of the standard estimator of the Pearson correlation coefficient, we applied a bias adjustment for $\hat{r}^2_{jk}$ using $\tilde{r}^2_{jk} = \frac{n-1}{n-2}\hat{r}^2_{jk} - \frac{1}{n-2}$.
* LD scores for variant $j$ were subsequently computed via $l_j = \sum_k \tilde{r}^2_{jk}$ with a radius of <b><u>1 MB</u></b>.
* For LDSC-compatible flat files, we only exported LD scores of high-quality HapMap 3 variants that are 1) in autosomes, 2) not in the MHC region, 3) biallelic SNPs, 4) with INFO > 0.9, and 5) MAF > 1% in UKB and gnomAD genome/exome (if available).
* We note that, since we applied covariate adjustment above, these LD scores are equivalent to the covariate-adjusted LD scores as described in [Luo, Y. & Li, X. et al., 2020](https://www.biorxiv.org/content/10.1101/503144v4)
