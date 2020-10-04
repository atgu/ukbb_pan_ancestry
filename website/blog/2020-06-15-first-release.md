---
id: first-release
title: First release
author: Pan UKBB Team
tags: [release, announcement]
---

We are thrilled to announce the release of GWAS summary statistics from the Pan-UK Biobank resource, which consists of genome-wide association analyses of 7,221 phenotypes across 6 continental ancestry groups in the [UK Biobank](https://www.ukbiobank.ac.uk/). Across all phenotype-ancestry pairs, we conducted 16,119 GWAS and meta-analyzed summary statistics for all available populations by trait. This release includes more than 20,000 individuals with primarily non-European ancestries, substantially increasing the diversity typically investigated in analyses of these data.

A summary of the breakdown included in this release is:

| Population   | Sample size | Total phenos | Categorical | Continuous | Phecode | ICD-10 | Biomarkers | Prescriptions |
|-------|-----------|----------------|------------|---------|-------|-------------|------------|---------------|
| AFR |      6636 |           2493 |        981 |     337 |   197 |        725  |         30 |           223 |
| AMR |       980 |           1105 |        423 |      31 |    20 |        561  |         30 |            40 |
| CSA |      8876 |           2771 |       1051 |     418 |   234 |        719  |         30 |           319 |
| EAS |      2709 |           1612 |        618 |      91 |    55 |        714  |         29 |           105 |
| EUR |    420531 |           7200 |       3672 |    1325 |   929 |        800  |         30 |           444 |
| MID |      1599 |           1372 |        509 |      83 |    52 |        591  |         30 |           107 |

<!--truncate-->

Rapidly developing this resource required a massively scalable computational framework, and we are thankful to the [Hail team](https://hail.is) for building and supporting the system that enabled this work. We developed a pipeline for these analyses using [Hail Batch](https://hail.is/docs/batch/index.html) and provide our [analytical code](https://github.com/atgu/ukbb_pan_ancestry) for reference to the community. We ran this pipeline in the Batch Service, a multi-tenant compute cluster in Google Cloud managed by the Hail team, which at the time of use, enabled the simultaneous use of up to 100,000 CPUs. Across all traits and phenotypes, the association tests required 3.8 million CPU-hours, which Hail Batch enabled to be completed in approximately 6 days (wall-clock).

We felt obligated to make this resource publicly available as soon as it was stable given its potential benefits to society, for example through contribution to pressing activities such as the [COVID-19 Host Genetics Initiative](https://www.covid19hg.org/). We will be further developing and analyzing this resource, with some additions forthcoming, as follows:
- We will be adding additional features to the website, including searchable and interactive summary statistics data as well as multi-ancestry visualizations of association results
- We will continually analyze and release updated COVID-19 phenotypes for all ancestry groups as they become available
- We are committed to computing LD scores in the UK Biobank for all ancestry groups analyzed in this resource and will provide these when they have been rigorously tested
- We will calculate heritabilities using LD score regression for each ancestry group and phenotype
- We will keep the phenotype manifest and release files updated with all further analyses

Diverse genetic studies are critical to equitably advance scientific discoveries and applications. However, history through current events teach us that some will unintentionally or willfully misuse this resource to advance racist agendas. Acknowledging the harm done to certain groups in the past in the name of science indicates the importance of careful communication of scientific research, its implications, and the intense vigilance required to ensure that disadvantaged groups are not further harmed by this and other related work. We have adopted several strategies in an effort to maximize benefits and minimize risks, most notably through the development of [a set of FAQs](https://pan.ukbb.broadinstitute.org/docs/background). These are designed to guide proper use and interpretation of these results.

Analyzing a more inclusive and diverse dataset increases power for discovery, enhances the resolution of these findings, and improves the generalizability of genetic associations across populations. We hope that this research will encourage future studies to make use of data that is traditionally left out of analysis and develop more diverse resources to ensure that genetics can execute on the mission of improving healthcare for all.
