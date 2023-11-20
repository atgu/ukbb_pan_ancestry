---
id: af-fix-sumstats-h2
title: Fixes to heritability estimates and adjustments to summary statisics schema
author: Rahul Gupta, on behalf of the Pan UKBB Team
tags: [sumstats, heritability, bugfix]
---

We have recently updated our summary statistic schema and our heritability estimates:

1. We have identified a bug in the computation of heritability point estimates involving use of improper allele frequencies. We have resolved the issue and accordingly recomputed all biobank-wide heritability analyses using all methods. These updated results can be found in the updated manifests.
2. With new heritability estimates, we have now recomputed QC for all summary statistics. This has resulted in a largely overlapping (but non-identical) set of 1091 QC-pass phenotype-ancestry pairs, with this new set used for QC-pass meta-analysis. Per-phenotype, per-ancestry association statistics remain identical.
3. We have recomputed the maximally independent set using the same approach as described previously, now incorporating the updated set of QC-pass phenotypes.
4. We have updated our summary statistics schema to now include clearly labeled -log10 p-values rather than ln p-values as previous. More information on the updated schema can be found on the [per-phenotype files](https://pan.ukbb.broadinstitute.org/docs/per-phenotype-files) page. Archived summary statistics using the previous schema can be found at  updated paths listed in the [archived sheet of the phenotype manifest](https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=268241601).

These changes should help improve the clarity and quality of our released data.