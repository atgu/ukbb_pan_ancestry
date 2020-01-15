import hail as hl

bucket = 'gs://ukb-diverse-pops'
REFERENCE_GENOME = 'GRCh37'
CHROMOSOMES = list(map(str, range(1, 23))) + ['X', 'XY']


def get_ukb_meta_pop_tsv_path():
    return f'{bucket}/pca/globalref_ukbb_pca_pops_rf_50.txt.bgz'


def get_ukb_meta():
    return hl.import_table(get_ukb_meta_pop_tsv_path(), key='s', impute=True)


def filter_pop(mt, pop):
    if pop != 'all': mt = mt.filter_cols(mt.pop == pop)
    return mt


def get_filtered_mt(pop: str = 'all'):
    mt = hl.read_matrix_table('gs://ukb31063/ukb31063.genotype.mt')
    meta_ht = get_ukb_meta()
    mt = mt.annotate_cols(**meta_ht.key_by(s=hl.str(meta_ht.s))[mt.s])
    return filter_pop(mt, pop)


def get_ukb_pheno_mt(pop: str = 'all'):
    from .phenotypes import get_ukb_pheno_mt_path
    path = get_ukb_pheno_mt_path()
    mt = hl.read_matrix_table(path)
    return filter_pop(mt, pop)

