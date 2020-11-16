---
id: ld-release
title: LD scores and matrices release
author: Masahiro Kanai, on behalf of the Pan UKBB Team
tags: [release, announcement, LD]
---

We are excited to announce [the release of genome-wide LD scores and matrices](https://pan.ukbb.broadinstitute.org/downloads) from the Pan-UK Biobank resource, which contains LD scores for (non-partitioned) [LD score regression](https://github.com/bulik/ldsc) analysis and full in-sample LD matrices for more extensive analyses such as fine-mapping.

<!--truncate-->

Computing large-scale LD matrices is always challenging: Thanks to the [Hail team](https://hail.is) for continuous development and support, we utilized Hail's block-distributed matrix implementation ([`BlockMatrix`](https://hail.is/docs/0.2/linalg/hail.linalg.BlockMatrix.html)) to enable massively scalable linear algebra on genotype matrix. With `BlockMatrix`, we computed genome-wide in-sample LD matrices (10 Mb radius) and scores for all six ancestry groups in approximately 16 hours (wall-clock) using 500 preemptible workers (n1-standard-8) on Google Cloud Dataproc (~64,000 CPU-hours).

Since our analysis included ancestry groups that are admixed and have more complex population structure, we need to account for that in downstream analyses such as LD score regression and fine-mapping. To this end, we applied covariate adjustment on genotype matrix for LD computation using the first 10 PCs as well as the other covariates used in the Pan-UKB GWAS ($age$, $sex$, $age*sex$, $age^2$, $age^2*sex$). Please find more details in [our technical documentation](https://pan.ukbb.broadinstitute.org/docs/ld).

Another potential concern might be using imputation-based LD scores for LD score regression, which is generally not recommended by the authors ([Bulik-Sullivan, BK. et al., 2015](https://www.nature.com/articles/ng.3211)). We rationalized our approach because not every ancestry group in our analysis has a decent sequence-based LD score available; but we also compared our UKBB LD scores (imputation-based) with gnomAD LD scores (sequence-based) for those populations available in gnomAD, and observed high concordance between the two LD scores. Please find our blog post [here](https://pan.ukbb.broadinstitute.org/blog/2020/10/29/ld-scores) by Rahul.

With these LD resources, we will continue to further analyze our GWAS results, including:
- LD score regression analysis to estimate heritabilities
- Fine-mapping analysis to identity causal variants of well-powered complex traits

We hope this additional release will encourage researchers to analyze diverse ancestry groups and/or develop statistical methods using population-matched in-sample LD.
