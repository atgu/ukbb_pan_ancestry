---
id: hail-format
title: Hail Format
---

## Release files

The results of this analysis are released in two main files on Google Cloud Storage (file format compatible with Hail >= 0.2.42):

- Summary statistics MatrixTable: `gs://ukb-diverse-pops-public/sumstats_release/results_full.mt` (12.78 T)
- Meta-analysis MatrixTable: `gs://ukb-diverse-pops-public/sumstats_release/meta_analysis.mt` (12.54 T)

These are also available on Amazon S3:

- Summary statistics MatrixTable: `s3://pan-ukb-us-east-1/sumstats_release/results_full.mt` (12.78 T)
- Meta-analysis MatrixTable: `s3://pan-ukb-us-east-1/sumstats_release/meta_analysis.mt` (12.54 T)

In addition, in-sample full LD matrices and scores are available on Amazon S3:
- LD BlockMatrix `s3://pan-ukb-us-east-1/ld_release/UKBB.{pop}.ldadj.bm` (43.3 T in total)
    - Size by population: AFR: 12.0 T, AMR: 3.3 T, CSA: 6.4T, EAS: 2.6T, EUR: 14.1T, MID: 4.9T
- Variant index Hail Table `s3://pan-ukb-us-east-1/ld_release/UKBB.{pop}.ldadj.variant.ht` (1.7 G in total)
- LD score Hail Table `s3://pan-ukb-us-east-1/ld_release/UKBB.{pop}.ldscore.ht` (4.0 G in total)

where `{pop}` represents one of the population abbreviations (i.e., AFR, AMR, CSA, EAS, EUR, or MID).

### Requester pays

Note that the files in the Google Cloud Storage bucket are "requester pays." In order to compute over these files or download them, you will need to specify a project which may be billed for access and download costs. The data are stored in a US multi-region bucket: thus, access to the dataset is free for use for Compute Engine instances started within US regions, as well as for full downloads within the US and Canada. When performing large analyses on the dataset, we suggest "bringing the compute to the data" and starting a VM or Dataproc cluster in a US region. You can browse the directory structure in a requester pays bucket with the `-u` flag (and note the `hl.init` call below to access the data using Hail):

```
gsutil -u your_project_id ls gs://ukb-diverse-pops-public/sumstats_release
```

## Using the libraries and files

The files on Google Cloud Platform can be accessed by cloning the [ukbb_pan_ancestry](https://github.com/atgu/ukbb_pan_ancestry) and the [ukb_common](https://github.com/Nealelab/ukb_common) repos and accessing them programmatically. We recommend using these functions, as they apply our QC metrics (e.g. the raw file contains 7,271 phenotypes, but use of this function will return 7,221 phenotypes after removing low-quality ones) and include convenience metrics such as lambda GC.

```
%%bash
git clone https://github.com/atgu/ukbb_pan_ancestry
git clone https://github.com/Nealelab/ukb_common
```


```
from ukbb_pan_ancestry import *

hl.init(spark_conf={'spark.hadoop.fs.gs.requester.pays.mode': 'AUTO',
                    'spark.hadoop.fs.gs.requester.pays.project.id': 'your_project_id'})

mt = load_final_sumstats_mt()

mt.describe()
```

## Results schema
The basic summary statistics have the following schema:

```
----------------------------------------
Column fields:
    'trait_type': str
    'phenocode': str
    'pheno_sex': str
    'coding': str
    'modifier': str
    'pheno_data': struct {
        n_cases: int32,
        n_controls: int32,
        heritability: float64,
        saige_version: str,
        inv_normalized: bool,
        pop: str,
        lambda_gc: float64,
        n_variants: int64,
        n_sig_variants: int64
    }
    'description': str
    'description_more': str
    'coding_description': str
    'category': str
    'n_cases_full_cohort_both_sexes': int64
    'n_cases_full_cohort_females': int64
    'n_cases_full_cohort_males': int64
----------------------------------------
Row fields:
    'locus': locus<GRCh37>
    'alleles': array<str>
    'rsid': str
    'varid': str
    'vep': struct {
        ...
    }
    'freq': array<struct {
        pop: str,
        ac: float64,
        af: float64,
        an: int64,
        gnomad_exomes_ac: int32,
        gnomad_exomes_af: float64,
        gnomad_exomes_an: int32,
        gnomad_genomes_ac: int32,
        gnomad_genomes_af: float64,
        gnomad_genomes_an: int32
    }>
    'pass_gnomad_exomes': bool
    'pass_gnomad_genomes': bool
    'n_passing_populations': int32
    'high_quality': bool
    'nearest_genes': array<struct {
        gene_id: str,
        gene_name: str,
        within_gene: bool
    }>
    'info': float64
----------------------------------------
Entry fields:
    'summary_stats': struct {
        AF_Allele2: float64,
        imputationInfo: float64,
        BETA: float64,
        SE: float64,
        `p.value.NA`: float64,
        `AF.Cases`: float64,
        `AF.Controls`: float64,
        Pvalue: float64,
        low_confidence: bool
    }
----------------------------------------
Column key: ['trait_type', 'phenocode', 'pheno_sex', 'coding', 'modifier']
Row key: ['locus', 'alleles']
----------------------------------------
```

### Columns (phenotypes)

The columns are indexed by phenotype using a composite key of trait type, phenocode, pheno_sex, coding, and modifier.  Trait types have one of the values below. `phenocode` typically corresponds to the Field from UK Biobank, or the specific ICD code or phecode, or a custom moniker. `pheno_sex` designates which sexes were run, and is marked as `both_sexes` for most traits, though some phecodes were restricted to `females` or `males`. The `coding` field is primarily used for categorical variables, to indicate which one-hot encoding was used (e.g. [coding 2](http://biobank.ctsu.ox.ac.uk/showcase/coding.cgi?id=100434) for [field 1747](http://biobank.ctsu.ox.ac.uk/showcase/field.cgi?id=1747)). Finally, `modifier` refers to any downstream modifications of the phenotype (e.g. `irnt` for inverse-rank normal transformation).

By default, the MatrixTable loaded by `load_final_sumstats_mt` returns one column per phenotype-population pair. We can see the number of unique phenotypes for each `trait_type` by:
```
phenotype_ht = mt.cols().collect_by_key()  # Converting into one element per phenotype
phenotype_ht.group_by('trait_type').aggregate(n_phenos=hl.agg.count()).show()
+-----------------+----------+
| trait_type      | n_phenos |
+-----------------+----------+
| "biomarkers"    |       30 |
| "categorical"   |     3684 |
| "continuous"    |      820 |
| "icd10"         |      915 |
| "phecode"       |     1327 |
| "prescriptions" |      445 |
+-----------------+----------+
```
You can explore the population-level data in more detail using:
```
phenotype_ht = mt.cols()
phenotype_ht.show(truncate=40, width=85)
+--------------+-----------+--------------+--------+----------+--------------------+
| trait_type   | phenocode | pheno_sex    | coding | modifier | pheno_data.n_cases |
+--------------+-----------+--------------+--------+----------+--------------------+
| str          | str       | str          | str    | str      |              int32 |
+--------------+-----------+--------------+--------+----------+--------------------+
| "biomarkers" | "30600"   | "both_sexes" | ""     | "irnt"   |               5759 |
| "biomarkers" | "30600"   | "both_sexes" | ""     | "irnt"   |                856 |
| "biomarkers" | "30600"   | "both_sexes" | ""     | "irnt"   |               7694 |
| "biomarkers" | "30600"   | "both_sexes" | ""     | "irnt"   |               2340 |
| "biomarkers" | "30600"   | "both_sexes" | ""     | "irnt"   |             367192 |
| "biomarkers" | "30600"   | "both_sexes" | ""     | "irnt"   |               1364 |
| "biomarkers" | "30610"   | "both_sexes" | ""     | "irnt"   |               6216 |
| "biomarkers" | "30610"   | "both_sexes" | ""     | "irnt"   |                938 |
| "biomarkers" | "30610"   | "both_sexes" | ""     | "irnt"   |               8422 |
| "biomarkers" | "30610"   | "both_sexes" | ""     | "irnt"   |               2572 |
+--------------+-----------+--------------+--------+----------+--------------------+
+-----------------------+-------------------------+--------------------------+
| pheno_data.n_controls | pheno_data.heritability | pheno_data.saige_version |
+-----------------------+-------------------------+--------------------------+
|                 int32 |                 float64 | str                      |
+-----------------------+-------------------------+--------------------------+
|                    NA |                2.54e-01 | "SAIGE_0.36.4"           |
|                    NA |                1.13e-01 | "SAIGE_0.36.4"           |
|                    NA |                2.41e-01 | "SAIGE_0.36.4"           |
|                    NA |                6.13e-02 | "SAIGE_0.36.4"           |
|                    NA |                6.45e-02 | "SAIGE_0.36.4"           |
|                    NA |                2.05e-01 | "SAIGE_0.36.4"           |
|                    NA |                2.85e-01 | "SAIGE_0.36.4"           |
|                    NA |                2.02e-01 | "SAIGE_0.36.4"           |
|                    NA |                3.92e-01 | "SAIGE_0.36.4"           |
|                    NA |                0.00e+00 | "SAIGE_0.36.4"           |
+-----------------------+-------------------------+--------------------------+
+---------------------------+----------------+----------------------+
| pheno_data.inv_normalized | pheno_data.pop | pheno_data.lambda_gc |
+---------------------------+----------------+----------------------+
|                      bool | str            |              float64 |
+---------------------------+----------------+----------------------+
|                     false | "AFR"          |             1.03e+00 |
|                     false | "AMR"          |             9.95e-01 |
|                     false | "CSA"          |             1.01e+00 |
|                     false | "EAS"          |             9.85e-01 |
|                     false | "EUR"          |             1.39e+00 |
|                     false | "MID"          |             1.00e+00 |
|                     false | "AFR"          |             1.03e+00 |
|                     false | "AMR"          |             9.98e-01 |
|                     false | "CSA"          |             1.07e+00 |
|                     false | "EAS"          |             9.67e-01 |
+---------------------------+----------------+----------------------+
+-----------------------+---------------------------+------------------------+
| pheno_data.n_variants | pheno_data.n_sig_variants | description            |
+-----------------------+---------------------------+------------------------+
|                 int64 |                     int64 | str                    |
+-----------------------+---------------------------+------------------------+
|              19174891 |                         1 | "Albumin"              |
|               9660572 |                         0 | "Albumin"              |
|              12514012 |                         7 | "Albumin"              |
|               8702282 |                         0 | "Albumin"              |
|              21020499 |                     38622 | "Albumin"              |
|              12296019 |                         0 | "Albumin"              |
|              19351952 |                       399 | "Alkaline phosphatase" |
|               9905980 |                         3 | "Alkaline phosphatase" |
|              12681888 |                       783 | "Alkaline phosphatase" |
|               8859024 |                       131 | "Alkaline phosphatase" |
+-----------------------+---------------------------+------------------------+
+------------------+--------------------+------------------------------------------+
| description_more | coding_description | category                                 |
+------------------+--------------------+------------------------------------------+
| str              | str                | str                                      |
+------------------+--------------------+------------------------------------------+
| NA               | NA                 | "Biological samples > Assay results >... |
| NA               | NA                 | "Biological samples > Assay results >... |
| NA               | NA                 | "Biological samples > Assay results >... |
| NA               | NA                 | "Biological samples > Assay results >... |
| NA               | NA                 | "Biological samples > Assay results >... |
| NA               | NA                 | "Biological samples > Assay results >... |
| NA               | NA                 | "Biological samples > Assay results >... |
| NA               | NA                 | "Biological samples > Assay results >... |
| NA               | NA                 | "Biological samples > Assay results >... |
| NA               | NA                 | "Biological samples > Assay results >... |
+------------------+--------------------+------------------------------------------+
+--------------------------------+-----------------------------+
| n_cases_full_cohort_both_sexes | n_cases_full_cohort_females |
+--------------------------------+-----------------------------+
|                          int64 |                       int64 |
+--------------------------------+-----------------------------+
|                         422605 |                      208336 |
|                         422605 |                      208336 |
|                         422605 |                      208336 |
|                         422605 |                      208336 |
|                         422605 |                      208336 |
|                         422605 |                      208336 |
|                         461525 |                      229156 |
|                         461525 |                      229156 |
|                         461525 |                      229156 |
|                         461525 |                      229156 |
+--------------------------------+-----------------------------+
+---------------------------+
| n_cases_full_cohort_males |
+---------------------------+
|                     int64 |
+---------------------------+
|                    179998 |
|                    179998 |
|                    179998 |
|                    179998 |
|                    179998 |
|                    179998 |
|                    194910 |
|                    194910 |
|                    194910 |
|                    194910 |
+---------------------------+
showing top 10 rows
```

More information about the GWAS run is found in the `pheno_data` struct. By default, when loading using `load_final_sumstats_mt`, the best practice QC parameters are used, which removes traits with a lambda GC < 0.5 or > 2. If this is undesirable, use `load_final_sumstats_mt(filter_phenos=False)`.

### Rows (variants)

The rows are indexed by locus and alleles. Direct annotations can be found in the `vep` schema, but we also provide a `nearest_genes` annotation for ease of analysis. Additionally, variant QC annotations are provided in the `high_quality` field (which is filtered to by default using `load_final_sumstats_mt` and can be switched off in the `filter_variants` parameter in that function).

### Entries (association tests)

The entry fields house the summary statistics themselves. Note that there is a `low_confidence` annotation that indicates a possible low-quality association test (allele count in cases or controls <= 3, or overall minor allele count < 20).

The resulting dataset can be filtered and annotated as a standard Hail MatrixTable:
```
mt = mt.filter_cols((mt.trait_type == 'phecode') & (mt.lambda_gc > 0.9) & (mt.lambda_gc < 1.1))
```

## Meta-analysis files

The meta-analysis results are in a similarly structured file:
```
meta_mt = hl.read_matrix_table(get_meta_analysis_results_path())
```
Here, the results are provided in an array, which includes the all-available-population meta-analysis in the 0th element `meta_mt.meta_analysis[0]` and leave-one-out meta-analyses in the remainder of the array.
```
Entry fields:
    'meta_analysis': array<struct {
        BETA: float64,
        SE: float64,
        Pvalue: float64,
        Q: float64,
        Pvalue_het: float64,
        N: int32,
        N_pops: int32,
        AF_Allele2: float64,
        AF_Cases: float64,
        AF_Controls: float64
    }>
```

### Combining the datasets
We also provide a function to annotate the overall sumstats MatrixTable with the largest meta-analysis for that phenotype.
```
mt = load_final_sumstats_mt()
mt = annotate_mt_with_largest_meta_analysis(mt)
```
If your analysis requires the simultaneous analysis of summary statistics from multiple populations (and not the meta-analysis), you can load the data with a similar structure to the meta-analysis MatrixTable (one column per phenotype, with population information packed into an array of entries and columns) using `load_final_sumstats_mt(separate_columns_by_pop=False)`.

## LD matrices

The LD matrices are in [`BlockMatrix`](https://hail.is/docs/0.2/linalg/hail.linalg.BlockMatrix.html) format. Please refer to [Hail's documentation](https://hail.is/docs/0.2/linalg/hail.linalg.BlockMatrix.html) for available operations on `BlockMatrix`.
```
from hail.linalg import BlockMatrix
bm = BlockMatrix.read(get_ld_matrix_path(pop='AFR'))
```
We note that the LD matrices were sparsified to a upper triangle (all elements of the lower triangle were zeroed out using [`BlockMatrix.sparsify_triangle`](https://hail.is/docs/0.2/linalg/hail.linalg.BlockMatrix.html#hail.linalg.BlockMatrix.sparsify_triangle)).

### Variant indices

To determine which row/column corresponds to which variant, we provide variant indices for `BlockMatrix` in Hail `Table` format.
```
ht_idx = hl.read_table(get_ld_variant_index_path(pop='AFR'))
```

The variant indices table has the following schema and `idx` corresponds to a row/column index in `BlockMatrix`.
```
----------------------------------------
Global fields:
    'n_samples': int32
    'pop': str
----------------------------------------
Row fields:
    'locus': locus<GRCh37>
    'alleles': array<str>
    'idx': int64
----------------------------------------
Key: ['locus', 'alleles']
----------------------------------------
```

### Extracting a subset of LD matrix

To extract a subset of LD matrix, you first need to identify indices of your variants of interest. Here, we provide two examples:

```
# filter by interval
interval = hl.parse_locus_interval('1:51572000-52857000')
ht_idx = ht_idx.filter(interval.contains(ht_idx.locus))

# or filter by a list of variant IDs (e.g., 1:51572412:A:G)
ht = hl.import_table('/path/to/your/list')
ht = ht.transmute(**hl.parse_variant(ht.variant)).key_by('locus', 'alleles')
ht_idx = ht_idx.join(ht, 'inner')
```

Then, you can filter the LD matrix into a subset using [`BlockMatrix.filter`](https://hail.is/docs/0.2/linalg/hail.linalg.BlockMatrix.html#hail.linalg.BlockMatrix.filter):
```
idx = ht_idx.idx.collect()
bm = bm.filter(idx, idx)
```

### Exporting a LD matrix to a flat file

Finally, to export a LD matrix to a flat file (txt file), you can use [`BlockMatrix.export`](https://hail.is/docs/0.2/linalg/hail.linalg.BlockMatrix.html#hail.linalg.BlockMatrix.export):

```
# Note: when you apply any operation on BlockMatrix,
# you need to write it first to storage before export
bm.write('/path/to/tmp/bm', force_row_major=True)
BlockMatrix.export(
    '/path/to/tmp/bm',
    '/path/to/flat_file.bgz',
    delimiter=' '
)
```

If your matrix is small enough to fit on memory, you can also directly export it to numpy via [`BlockMatrix.to_numpy`](https://hail.is/docs/0.2/linalg/hail.linalg.BlockMatrix.html#hail.linalg.BlockMatrix.to_numpy).

```
np_mat = bm.to_numpy()
```

## LD scores

The LD scores are in Hail [Table](https://hail.is/docs/0.2/hail.Table.html) format. For LDSC-compatible flat files, you can find them [here](https://pan-ukb-us-east-1.s3.amazonaws.com/ld_release/UKBB.ALL.ldscore.tar.gz).
```
ht = hl.read_table(get_ld_score_ht_path(pop='AFR'))
```

The LD score table has the following schema.
```
----------------------------------------
Global fields:
    None
----------------------------------------
Row fields:
    'locus': locus<GRCh37>
    'alleles': array<str>
    'rsid': str
    'varid': str
    'AF': float64
    'ld_score': float64
----------------------------------------
Key: ['locus', 'alleles']
----------------------------------------
```
