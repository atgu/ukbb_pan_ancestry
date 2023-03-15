---
id: changelog
title: Changelog
---

### Version 0.4

- Bugfix implemented for incorrect allele frequencies used for heritability estimates, with all estimates now recomputed and all downstream QC and analyses updated.
- Updated sumstats schema now using $-\log_{10} p$-values with field names reflective of this change.

(released Mar 16th 2023)

### Version 0.3

- Added h2 estimates, summary statistics QC flags, and maximally independent set of QC-pass phenotypes. 
- Updated summary statistics format to report $\log p$-values, fixed traits showing standard errors of 0, and reproduced cross-ancestry meta-analyses to incorporate these updated results. 
- Produced new cross-ancestry meta-analyses using only QC-pass ancestry-trait pairs. 

(released Mar 17th 2022)

### Version 0.2
Added LD scores and matrices (released Oct 29, 2020).

### Version 0.1

Initial version (released June 16, 2020).
