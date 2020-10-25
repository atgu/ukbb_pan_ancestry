import hail as hl

bucket = 'gs://ukb-diverse-pops'
temp_bucket = 'gs://ukbb-diverse-temp-30day'
temp_bucket_7day = 'gs://ukbb-diverse-temp-7day'
public_bucket = 'gs://ukb-diverse-pops-public'
public_bucket_free = 'gs://ukb-diverse-pops-public-free'

REFERENCE_GENOME = 'GRCh37'
CHROMOSOMES = list(map(str, range(1, 23))) + ['X', 'XY']
POPS = ['AFR', 'AMR', 'CSA', 'EAS', 'EUR', 'MID']

MIN_CASES = 50
MIN_CASES_ALL = 100
MIN_CASES_EUR = 100


def get_hq_samples():
    ht = hl.import_table(f'{bucket}/misc/ukb31063_samples_qc_FULL.txt', no_header=True)
    drop_samples = hl.import_table(f'{bucket}/misc/ukb31063.withdrawn_samples_20190321.txt', no_header=True, key='f0')
    ht = ht.key_by(s=ht.f0).drop('f0')
    return ht.filter(hl.is_missing(drop_samples[ht.s]))


def get_pruned_tsv_path():
    return f'{bucket}/pca/ukb_diverse_pops_pruned.tsv.bgz'


def get_age_sex_tsv_path():
    return f'{bucket}/Phenotypes/uk_round2_allSamples_phenos_phesant.6148_5.tsv.gz'


def get_covariates_ht_path(extension: str = 'ht'):
    return f'{bucket}/pca/all_pops_non_eur_pruned_within_pop_pc_covs.{extension}'


def get_covariates(key_type = hl.str):
    ht = hl.read_table(get_covariates_ht_path())
    return ht.key_by(s=key_type(ht.s))


def get_final_sample_set():
    return f'{bucket}/misc/final_samples.txt.bgz'


def get_ukb_meta_pop_tsv_path():
    # WARNING: deprecated (original file for first iteration), remains for historical record
    return f'{bucket}/pca/globalref_ukbb_pca_pops_rf_50.txt.bgz'


def get_ukb_meta(key_type = hl.tstr):
    # WARNING: deprecated (original file for first iteration), remains for historical record
    return hl.import_table(get_ukb_meta_pop_tsv_path(), key='s', impute=True, types={'s': key_type})


def get_ukb_pheno_mt(pop: str = 'all'):
    from .phenotypes import get_ukb_pheno_mt_path
    mt = hl.read_matrix_table(get_ukb_pheno_mt_path())
    # mt = mt.annotate_rows(**get_ukb_meta(key_type=hl.tint32)[mt.row_key])
    mt = mt.annotate_rows(**get_covariates(key_type=hl.int32)[mt.row_key])
    if pop != 'all':
        mt = mt.filter_rows(mt.pop == pop)
    return mt
