import hail as hl

bucket = 'gs://ukb-diverse-pops'
REFERENCE_GENOME = 'GRCh37'
CHROMOSOMES = list(map(str, range(1, 23))) + ['X', 'XY']


def get_ukb_meta_pop_tsv_path():
    return f'{bucket}/globalref_ukbb_pca_pops_rf_50.txt.bgz'

ukb_imputed_bgen_path = 'gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183/imputed/ukb_imp_chr{}_v3.bgen'
ukb_imputed_info_path = 'gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183/imputed/ukb_mfi_chr{}_v3.txt'
ukb_imputed_info_ht_path = f'{bucket}/imputed/ukb_mfi_v3.ht'


def get_sample_file(chromosome: str = '1'):
    if chromosome not in ('X', 'XY'):
        chromosome = 'autosomes'
    return f'gs://ukb31063/ukb31063.{chromosome}.sample'


def get_ukb_imputed_data(chromosome: str = '1', variant_list: hl.Table = None, entry_fields = ('dosage', )):
    if chromosome == 'all':
        chromosome = '{' + ','.join(map(str, range(1, 23))) + '}'
    add_args = {}
    if variant_list is not None:
        add_args['variants'] = variant_list
    return hl.import_bgen(ukb_imputed_bgen_path.format(chromosome), entry_fields=entry_fields,
                          sample_file=get_sample_file(chromosome), **add_args)

ukb_af_ht_path = f'{bucket}/imputed/ukb_frequencies.ht'
