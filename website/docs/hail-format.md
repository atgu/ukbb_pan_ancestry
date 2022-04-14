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

The files on Google Cloud Platform can be accessed by cloning the [ukbb_pan_ancestry](https://github.com/atgu/ukbb_pan_ancestry) and the [ukb_common](https://github.com/Nealelab/ukb_common) repos and accessing them programmatically. We recommend using these functions, as they allow for automatic application of our QC metrics as well as inclusion of all [QC flags](https://pan.ukbb.broadinstitute.org/docs/qc#quality-control-of-summary-statistics) and convenience metrics such as lambda GC. By default, when loading using `load_final_sumstats_mt`, the best practice QC parameters are used, which removes traits with a lambda GC < 0.5 or > 5 as well as applying all [QC filters](https://pan.ukbb.broadinstitute.org/docs/qc#quality-control-of-summary-statistics). This results in importing summary statistics for 527 traits; if it is preferable to load all traits with exported summary statistics (e.g., only applying the lambda GC < 0.5 or > 5 filter), use `load_final_sumstats_mt(filter_pheno_h2_qc=False)`, resulting in 7,228 traits. If any filtering is undesirable, use `load_final_sumstats_mt(filter_pheno_h2_qc=False, filter_phenos=False)`, which will import all 7,271 traits.

```
%%bash
git clone https://github.com/atgu/ukbb_pan_ancestry
git clone https://github.com/Nealelab/ukb_common
```


```
from ukbb_pan_ancestry import *

hl.init(spark_conf={'spark.hadoop.fs.gs.requester.pays.mode': 'AUTO',
                    'spark.hadoop.fs.gs.requester.pays.project.id': 'your_project_id'})

# loads all results for which sumstats were exported (lambda GC < 0.5 or > 5)
mt = load_final_sumstats_mt(filter_pheno_h2_qc=False)
# use filter_pheno_h2_qc=True to filter to just ancestry-trait pairs passing all QC
# mt = load_final_sumstats_mt(filter_pheno_h2_qc=True)

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
        heritability: struct {
            estimates: struct {
                ldsc: struct {
                    h2_liability: float64,
                    h2_liability_se: float64,
                    h2_z: float64,
                    h2_observed: float64,
                    h2_observed_se: float64,
                    intercept: float64,
                    intercept_se: float64,
                    ratio: float64,
                    ratio_se: float64
                },
                sldsc_25bin: struct {
                    h2_liability: float64,
                    h2_liability_se: float64,
                    h2_z: float64,
                    h2_observed: float64,
                    h2_observed_se: float64,
                    intercept: float64,
                    intercept_se: float64,
                    ratio: float64,
                    ratio_se: float64
                },
                rhemc_25bin: struct {
                    h2_liability: float64,
                    h2_liability_se: float64,
                    h2_z: float64,
                    h2_observed: float64,
                    h2_observed_se: float64
                },
                rhemc_8bin: struct {
                    h2_liability: float64,
                    h2_liability_se: float64,
                    h2_observed: float64,
                    h2_observed_se: float64,
                    h2_z: float64
                },
                rhemc_25bin_50rv: struct {
                    h2_observed: float64,
                    h2_observed_se: float64,
                    h2_liability: float64,
                    h2_liability_se: float64,
                    h2_z: float64
                },
                final: struct {
                    h2_observed: float64,
                    h2_observed_se: float64,
                    h2_liability: float64,
                    h2_liability_se: float64,
                    h2_z: float64
                }
            },
            qcflags: struct {
                GWAS_run: bool,
                ancestry_reasonable_n: bool,
                defined_h2: bool,
                significant_z: bool,
                in_bounds_h2: bool,
                normal_lambda: bool,
                normal_ratio: bool,
                EUR_plus_1: bool,
                pass_all: bool
            },
            N_ancestry_QC_pass: int32
        },
        saige_version: str,
        inv_normalized: bool,
        pop: str,
        lambda_gc: float64,
        n_variants: int64,
        n_sig_variants: int64,
        saige_heritability: float64
    }
    'description': str
    'description_more': str
    'coding_description': str
    'category': str
    'n_cases_full_cohort_both_sexes': int64
    'n_cases_full_cohort_females': int64
    'n_cases_full_cohort_males': int64
    'pop_index': int32
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

# results for all exported sumstats
+-----------------+----------+
| trait_type      | n_phenos |
+-----------------+----------+
| str             |    int64 |
+-----------------+----------+
| "biomarkers"    |       30 |
| "categorical"   |     3686 |
| "continuous"    |      820 |
| "icd10"         |      921 |
| "phecode"       |     1326 |
| "prescriptions" |      445 |
+-----------------+----------+

# results for full QC pass-only sumstats
+-----------------+----------+
| trait_type      | n_phenos |
+-----------------+----------+
| str             |    int64 |
+-----------------+----------+
| "biomarkers"    |       23 |
| "categorical"   |      179 |
| "continuous"    |      206 |
| "icd10"         |       34 |
| "phecode"       |       64 |
| "prescriptions" |       21 |
+-----------------+----------+
```
You can explore the population-level data in more detail using (several fields removed for brevity):
```
phenotype_ht = mt.cols()
phenotype_ht.show(truncate=40, width=105)
+--------------+-----------+--------------+--------+----------+--------------------+
| trait_type   | phenocode | pheno_sex    | coding | modifier | pheno_data.n_cases |
+--------------+-----------+--------------+--------+----------+--------------------+
| str          | str       | str          | str    | str      |              int32 |
+--------------+-----------+--------------+--------+----------+--------------------+
| "biomarkers" | "30600"   | "both_sexes" | ""     | "irnt"   |               7694 |
| "biomarkers" | "30600"   | "both_sexes" | ""     | "irnt"   |             367192 |
| "biomarkers" | "30610"   | "both_sexes" | ""     | "irnt"   |               8422 |
| "biomarkers" | "30610"   | "both_sexes" | ""     | "irnt"   |             400988 |
| "biomarkers" | "30620"   | "both_sexes" | ""     | "irnt"   |               6214 |
| "biomarkers" | "30620"   | "both_sexes" | ""     | "irnt"   |               8407 |
| "biomarkers" | "30620"   | "both_sexes" | ""     | "irnt"   |             400822 |
| "biomarkers" | "30620"   | "both_sexes" | ""     | "irnt"   |               1499 |
| "biomarkers" | "30630"   | "both_sexes" | ""     | "irnt"   |               7679 |
| "biomarkers" | "30630"   | "both_sexes" | ""     | "irnt"   |             364987 |
+--------------+-----------+--------------+--------+----------+--------------------+

+-----------------------+------------------------------------------+
| pheno_data.n_controls | pheno_data.heritability.estimates.lds... |
+-----------------------+------------------------------------------+
|                 int32 |                                  float64 |
+-----------------------+------------------------------------------+
|                    NA |                                 1.62e-01 |
|                    NA |                                 1.18e-01 |
|                    NA |                                 1.98e-01 |
|                    NA |                                 2.18e-01 |
|                    NA |                                 1.27e-01 |
|                    NA |                                 1.48e-02 |
|                    NA |                                 1.14e-01 |
|                    NA |                                -2.67e-01 |
|                    NA |                                 1.30e-01 |
|                    NA |                                 1.89e-01 |
+-----------------------+------------------------------------------+

+------------------------------------------+------------------------------------------+
| pheno_data.heritability.qcflags.norma... | pheno_data.heritability.qcflags.norma... |
+------------------------------------------+------------------------------------------+
|                                     bool |                                     bool |
+------------------------------------------+------------------------------------------+
|                                     True |                                     True |
|                                     True |                                     True |
|                                     True |                                     True |
|                                     True |                                     True |
|                                     True |                                     True |
|                                     True |                                     True |
|                                     True |                                     True |
|                                     True |                                     True |
|                                     True |                                     True |
|                                     True |                                     True |
+------------------------------------------+------------------------------------------+

+------------------------------------------+------------------------------------------+
| pheno_data.heritability.qcflags.EUR_p... | pheno_data.heritability.qcflags.pass_all |
+------------------------------------------+------------------------------------------+
|                                     bool |                                     bool |
+------------------------------------------+------------------------------------------+
|                                     True |                                     True |
|                                     True |                                     True |
|                                     True |                                     True |
|                                     True |                                     True |
|                                     True |                                     True |
|                                     True |                                     True |
|                                     True |                                     True |
|                                     True |                                     True |
|                                     True |                                     True |
|                                     True |                                     True |
+------------------------------------------+------------------------------------------+

+------------------------------------------+--------------------------+---------------------------+
| pheno_data.heritability.N_ancestry_QC... | pheno_data.saige_version | pheno_data.inv_normalized |
+------------------------------------------+--------------------------+---------------------------+
|                                    int32 | str                      |                      bool |
+------------------------------------------+--------------------------+---------------------------+
|                                        2 | "SAIGE_0.36.4"           |                     False |
|                                        2 | "SAIGE_0.36.4"           |                     False |
|                                        2 | "SAIGE_0.36.4"           |                     False |
|                                        2 | "SAIGE_0.44.5"           |                     False |
|                                        4 | "SAIGE_0.36.4"           |                     False |
|                                        4 | "SAIGE_0.36.4"           |                     False |
|                                        4 | "SAIGE_0.44.5"           |                     False |
|                                        4 | "SAIGE_0.36.4"           |                     False |
|                                        3 | "SAIGE_0.36.4"           |                     False |
|                                        3 | "SAIGE_0.44.5"           |                     False |
+------------------------------------------+--------------------------+---------------------------+

+----------------+----------------------+-----------------------+---------------------------+
| pheno_data.pop | pheno_data.lambda_gc | pheno_data.n_variants | pheno_data.n_sig_variants |
+----------------+----------------------+-----------------------+---------------------------+
| str            |              float64 |                 int64 |                     int64 |
+----------------+----------------------+-----------------------+---------------------------+
| "CSA"          |             1.03e+00 |              12200078 |                         6 |
| "EUR"          |             1.37e+00 |              20561726 |                     37450 |
| "CSA"          |             1.02e+00 |              12364741 |                       772 |
| "EUR"          |             1.67e+00 |              20739978 |                     89683 |
| "AFR"          |             1.02e+00 |              18630599 |                        77 |
| "CSA"          |             1.02e+00 |              12362444 |                         0 |
| "EUR"          |             1.42e+00 |              20739238 |                     38220 |
| "MID"          |             9.89e-01 |              12328418 |                         0 |
| "CSA"          |             1.01e+00 |              12195348 |                       378 |
| "EUR"          |             1.63e+00 |              20547047 |                     62484 |
+----------------+----------------------+-----------------------+---------------------------+

+-------------------------------+----------------------------+------------------+--------------------+
| pheno_data.saige_heritability | description                | description_more | coding_description |
+-------------------------------+----------------------------+------------------+--------------------+
|                       float64 | str                        | str              | str                |
+-------------------------------+----------------------------+------------------+--------------------+
|                      2.41e-01 | "Albumin"                  | NA               | NA                 |
|                      6.45e-02 | "Albumin"                  | NA               | NA                 |
|                      3.92e-01 | "Alkaline phosphatase"     | NA               | NA                 |
|                      1.31e-01 | "Alkaline phosphatase"     | NA               | NA                 |
|                      3.94e-01 | "Alanine aminotransferase" | NA               | NA                 |
|                      2.19e-01 | "Alanine aminotransferase" | NA               | NA                 |
|                      6.28e-02 | "Alanine aminotransferase" | NA               | NA                 |
|                      2.05e-01 | "Alanine aminotransferase" | NA               | NA                 |
|                      3.58e-01 | "Apolipoprotein A"         | NA               | NA                 |
|                      1.22e-01 | "Apolipoprotein A"         | NA               | NA                 |
+-------------------------------+----------------------------+------------------+--------------------+

+------------------------------------------+--------------------------------+
| category                                 | n_cases_full_cohort_both_sexes |
+------------------------------------------+--------------------------------+
| str                                      |                          int64 |
+------------------------------------------+--------------------------------+
| "Biological samples > Assay results >... |                         422605 |
| "Biological samples > Assay results >... |                         422605 |
| "Biological samples > Assay results >... |                         461525 |
| "Biological samples > Assay results >... |                         461525 |
| "Biological samples > Assay results >... |                         461326 |
| "Biological samples > Assay results >... |                         461326 |
| "Biological samples > Assay results >... |                         461326 |
| "Biological samples > Assay results >... |                         461326 |
| "Biological samples > Assay results >... |                         420088 |
| "Biological samples > Assay results >... |                         420088 |
+------------------------------------------+--------------------------------+

+-----------------------------+---------------------------+-----------+
| n_cases_full_cohort_females | n_cases_full_cohort_males | pop_index |
+-----------------------------+---------------------------+-----------+
|                       int64 |                     int64 |     int32 |
+-----------------------------+---------------------------+-----------+
|                      208336 |                    179998 |         0 |
|                      208336 |                    179998 |         1 |
|                      229156 |                    194910 |         0 |
|                      229156 |                    194910 |         1 |
|                      229118 |                    194764 |         0 |
|                      229118 |                    194764 |         1 |
|                      229118 |                    194764 |         2 |
|                      229118 |                    194764 |         3 |
|                      206413 |                    179623 |         0 |
|                      206413 |                    179623 |         1 |
+-----------------------------+---------------------------+-----------+
showing top 10 rows
```

More information about the GWAS run is found in the `pheno_data` struct. This struct also includes all heritability information and QC flags. More details on heritability estimation methods are forthcoming, but can be previewed [here](https://pan.ukbb.broadinstitute.org/blog/2022/04/11/h2-qc-updated-sumstats). Descriptions of the heritability fields can be found [below](https://pan.ukbb.broadinstitute.org/docs/hail-format#heritability) and more information on QC flags can be found [here.](https://pan.ukbb.broadinstitute.org/docs/qc#quality-control-of-summary-statistics)

### Rows (variants)

The rows are indexed by locus and alleles. Direct annotations can be found in the `vep` schema, but we also provide a `nearest_genes` annotation for ease of analysis. Additionally, variant QC annotations are provided in the `high_quality` field (which is filtered to by default using `load_final_sumstats_mt` and can be switched off in the `filter_variants` parameter in that function).

### Entries (association tests)

:::note
Please note that p-values are now stored as log p-values to avoid underflow.
:::

The entry fields house the summary statistics themselves. Note that there is a `low_confidence` annotation that indicates a possible low-quality association test (allele count in cases or controls <= 3, or overall minor allele count < 20).

The resulting dataset can be filtered and annotated as a standard Hail MatrixTable:
```
mt = mt.filter_cols((mt.trait_type == 'phecode') & (mt.lambda_gc > 0.9) & (mt.lambda_gc < 1.1))
```

## Meta-analysis files

:::note
Please note that p-values are now stored as log p-values to avoid underflow.
:::

The meta-analysis results are in a similarly structured file which can be obtained as such:
```
meta_mt = hl.read_matrix_table(get_meta_analysis_results_path())
```
By default, this function imports the "high-quality" meta-analyses, which are meta-analyses of only ancestry groups passing all [QC filters](https://pan.ukbb.broadinstitute.org/docs/qc#quality-control-of-summary-statistics). Naturally these meta-analyses are only available for those phenotypes for which at least two ancestries passed all QC (since a requirement for QC PASS is that a trait passes QC in EUR and at least one other ancestry). If interested in meta-analyses of all ancestries for which GWAS was run for a given phenotype, import using `hl.read_matrix_table(get_meta_analysis_results_path(filter_pheno_h2_qc=False))`.

Both versions of the meta-analysis table have results provided in an array, which includes the all-available-population meta-analysis in the 0th element `meta_mt.meta_analysis[0]` and leave-one-out meta-analyses in the remainder of the array.
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

## Heritability estimates

The heritability estimates can be found as a [flat file manifest](https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=1797288938) or as a Hail [Table](https://hail.is/docs/0.2/hail.Table.html).
```
ht = hl.read_table(get_h2_ht())
```

This table has the following schema:
```
----------------------------------------
Global fields:
    None
----------------------------------------
Row fields:
    'trait_type': str
    'phenocode': str
    'pheno_sex': str
    'coding': str
    'modifier': str
    'heritability': array<struct {
        pop: str,
        estimates: struct {
            ldsc: struct {
                h2_liability: float64,
                h2_liability_se: float64,
                h2_z: float64,
                h2_observed: float64,
                h2_observed_se: float64,
                intercept: float64,
                intercept_se: float64,
                ratio: float64,
                ratio_se: float64
            },
            sldsc_25bin: struct {
                h2_liability: float64,
                h2_liability_se: float64,
                h2_z: float64,
                h2_observed: float64,
                h2_observed_se: float64,
                intercept: float64,
                intercept_se: float64,
                ratio: float64,
                ratio_se: float64
            },
            rhemc_25bin: struct {
                h2_liability: float64,
                h2_liability_se: float64,
                h2_z: float64,
                h2_observed: float64,
                h2_observed_se: float64
            },
            rhemc_8bin: struct {
                h2_liability: float64,
                h2_liability_se: float64,
                h2_observed: float64,
                h2_observed_se: float64,
                h2_z: float64
            },
            rhemc_25bin_50rv: struct {
                h2_observed: float64,
                h2_observed_se: float64,
                h2_liability: float64,
                h2_liability_se: float64,
                h2_z: float64
            },
            final: struct {
                h2_observed: float64,
                h2_observed_se: float64,
                h2_liability: float64,
                h2_liability_se: float64,
                h2_z: float64
            }
        },
        qcflags: struct {
            GWAS_run: bool,
            ancestry_reasonable_n: bool,
            defined_h2: bool,
            significant_z: bool,
            in_bounds_h2: bool,
            normal_lambda: bool,
            normal_ratio: bool,
            EUR_plus_1: bool,
            pass_all: bool
        },
        N_ancestry_QC_pass: int32
    }>
----------------------------------------
Key: ['trait_type', 'phenocode', 'pheno_sex', 'coding', 'modifier']
----------------------------------------
```

Note that this is very similar to the heritability struct in the [results schema](#results-schema) -- the `load_final_sumstats_mt()` automatically uses `get_h2_ht()` to import the heritability table and annotate it into the column schema of the summary statistics results table.

Here the heritability field is an `array` of structs, where the elements of the array represent the ancestries for which heritability estimates are available. The corresponding ancestry is found in the `heritability.pop` field. To obtain one row for each ancestry-trait pair, use the following commmand:
```
ht = ht.explode('heritability')
```

The `heritability.estimates` struct contains point estimates and significance test results for:

- Univariate LD score regression (`ldsc`), run using [LD score flat files](https://pan.ukbb.broadinstitute.org/docs/ld#ld-scores) using high-quality HapMap3 SNPs with MAF $\geq$ 0.01 with summary statistics exported from the [results table](#results-schema).
- Stratified LD score regression (`sldsc_25bin`), run using the same summary statistcs as `ldsc` with LD scores generated from SNPs in 5 MAF and 5 LD score bins.
- Randomized Haseman-Elston (`rhemc_8bin`) using genotype data with 2 MAF and 4 LD score bins using default settings. This was predominantly used to analyze non-EUR ancestry groups but includes a small set of estimates for traits in EUR.
- Randomized Haseman-Elston (`rhemc_25bin`) using genotype data with 5 MAF and 5 LD score bins using default settings. Only run for non-EUR ancestry groups.
- Randomized Haseman-Elston (`rhemc_25bin_50rv`) using genotype data with 5 MAF and 5 LD score bins using 50 random variables to reduce run-to-run variability at a slightly higher computational cost. Only run for non-EUR ancestry groups.

Within each of the above methods, we produce observed- and liability-scale estimates as well as standard errors and the z-score for the test of $h^2 > 0$. We also produce LDSC intercept and ratio estimates when relevant.

The `heritability.qcflags` struct contains the results from our sequential QC filtering scheme -- see [quality control](https://pan.ukbb.broadinstitute.org/docs/qc#quality-control-of-summary-statistics) for more details. `N_ancestry_QC_pass ` refers to the number of ancestries passing all QC for a given trait.

More information can be found [here](https://pan.ukbb.broadinstitute.org/docs/heritability) on the heritability estimation approach, and important caveats can be found [here](https://pan.ukbb.broadinstitute.org/docs/qc#heritability).