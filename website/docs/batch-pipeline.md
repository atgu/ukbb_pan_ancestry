---
id: batch-pipeline
title: Batch Pipeline
---

## Batch pipeline

### Initial run

- Phenotypes are loaded into a MatrixTable using [load_phenotype_data.py](https://github.com/atgu/ukbb_pan_ancestry/blob/master/load_phenotype_data.py)

- The pipeline is run using Hail Batch in the script [saige_pan_ancestry.py](https://github.com/atgu/ukbb_pan_ancestry/blob/master/saige_pan_ancestry.py)

- Results are loaded into a MatrixTable using [load_all_results.py](https://github.com/atgu/ukbb_pan_ancestry/blob/master/load_all_results.py)

For the initial run, we ran the Batch submission script ([saige_pan_ancestry.py](https://github.com/atgu/ukbb_pan_ancestry/blob/master/saige_pan_ancestry.py)) on a VM on the cloud as it requires lots of memory and fast access to `gs://` files. The initial pipeline spanned approximately 4 million CPU-hours.

### How to run additional phenotypes

*Internal use only*

Custom phenotypes are the easiest to run and add. A file with:
```
eid pheno1name pheno2name
1234 0 1
1235 1 1
```
uploaded to `gs://ukb-diverse-pops/Phenotypes/Everyone/custom/TRAITTYPE_NAME.txt` (where `TRAITTYPE` is one of `categorical` or `continuous`).

#### Step 1 - Load phenotype data

This file can then be loaded using:
```
python load_phenotype_data.py --add_dataset TRAITTYPE_NAME
```

If you do not have the google cloud connector installed, you can do this as a submission (`hailctl dataproc submit load_phenotype_data.py --add_dataset TRAITTYPE_NAME`); note that in either case, you will need to have the `ukbb_common` and `ukbb_pan_ancestry` packages on your path. This script will make a backup copy of the existing phenotype MatrixTable and overwrite the master file with a new one including the new phenotype(s).

#### Step 2 - Run tests

You can then run the phenotype through the pipeline using batch (can do this locally since jobs are dispatched to the cloud anyway).

```
python saige_pan_ancestry.py --run_all_phenos --specific_phenos '.*-NAME'
```

will run all the phenotypes matching that regex. You can add `--dry_run` in order to first test the pipeline to see how many phenotypes will be run (there will be lots of output but any line starting with `Submitting` will tell you how many tasks of each type will be submit - there should be one for each phenotype for most tasks, and ~600 (~3000 for EUR) per phenotype for `saige_task`). Minimum case count can be overrode with `--skip_case_count_filter`.

When the pipeline completes, the results will have been loaded into Hail Tables at 

```gs://ukb-diverse-pops/results/result/POP/TRAITTYPE-pheno1-NAME/variant_results.ht```

This directory will also have `.txt` and `.log` files for each chunk run (can safely delete the `.txt` files at this point). Note that not all phenotypes may complete (especially if a SAIGE null model fails to converge).

#### Step 3 - Load into Results MatrixTable

Once in a while, to regenerate the MatrixTable, run:

```
python load_all_results.py --run_additional_load --overwrite
```

As above, `--dry_run` will give you an indication of how many new phenotypes will be loaded (for AFR). As with the phenotype data, a backup will be made of the existing population-specific MatrixTables and then the full MatrixTable will be overwritten.
