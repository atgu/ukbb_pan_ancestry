import hail as hl
from .generic import *

ukb_imputed_bgen_path = 'gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183/imputed/ukb_imp_chr{}_v3.bgen'
ukb_imputed_info_path = 'gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183/imputed/ukb_mfi_chr{}_v3.txt'
ukb_imputed_info_ht_path = f'{bucket}/imputed/ukb_mfi_v3.ht'


def get_sample_file(chromosome: str = '1'):
    if chromosome not in ('X', 'XY'):
        chromosome = 'autosomes'
    elif not chromosome.startswith('chr'):
        chromosome = f'chr{chromosome}'
    return f'gs://ukb31063/ukb31063.{chromosome}.sample'


def get_ukb_imputed_data(chromosome: str = '1', variant_list: hl.Table = None, entry_fields = ('GP', )):
    if chromosome == 'all':
        chromosome = '{' + ','.join(map(str, range(1, 23))) + '}'
    add_args = {}
    if variant_list is not None:
        add_args['variants'] = variant_list
    return hl.import_bgen(ukb_imputed_bgen_path.format(chromosome), entry_fields=entry_fields,
                          sample_file=get_sample_file(chromosome), **add_args)


def get_filtered_mt(chrom: str = 'all',
                    pop: str = 'all',
                    imputed: bool = True,
                    min_mac: int = 20,
                    entry_fields=('GP',),
                    filter_mac_instead_of_ac: bool = False):

    # get ac or mac based on filter_mac_instead_of_ac
    def get_ac(af, an):
        if filter_mac_instead_of_ac:
            # Note that the underlying file behind get_ukb_af_ht_path() accidentally double af and halve an
            return (1.0 - hl.abs(1.0 - af)) * an
        else:
            return af * an

    if imputed:
        ht = hl.read_table(get_ukb_af_ht_path())
        if pop == 'all':
            ht = ht.filter(hl.any(lambda x: get_ac(ht.af[x], ht.an[x]) >= min_mac, hl.literal(POPS)))
        else:
            ht = ht.filter(get_ac(ht.af[pop], ht.an[pop]) >= min_mac)
        mt = get_ukb_imputed_data(chrom, variant_list=ht, entry_fields=entry_fields)
    else:
        mt = hl.read_matrix_table('gs://ukb31063/ukb31063.genotype.mt')

    covariates_ht = get_covariates()
    hq_samples_ht = get_hq_samples()
    mt = mt.annotate_cols(**covariates_ht[mt.s])
    mt = mt.filter_cols(hl.is_defined(mt.pop) & hl.is_defined(hq_samples_ht[mt.s]))

    if pop != 'all': mt = mt.filter_cols(mt.pop == pop)
    return mt


def get_filtered_mt_with_x(pop: str = 'all'):
    mt = get_filtered_mt(pop=pop, entry_fields=('dosage',))
    mt_x = get_filtered_mt('X', pop=pop, entry_fields=('dosage',))
    # assert mt.s.collect() == mt_x.s.collect()
    mt = mt.union_rows(mt_x)
    return mt


def get_ukb_af_ht_path(with_x = True, repart=False):
    return f'{bucket}/imputed/ukb_frequencies{"_with_x" if with_x else ""}{".repart" if repart else ""}.ht'


def get_ukb_vep_path():
    return f'{bucket}/results/misc/ukb.vep.ht'


def get_ukb_grm_mt_path(pop: str, data_iteration: int = 0):
    suffix = f'.{data_iteration}' if data_iteration else ""
    return f'{bucket}/results/misc/ukb.{pop}.for_grm{suffix}.mt'


def get_ukb_grm_pruned_ht_path(pop: str, window_size: str = '1e6'):
    cut = '' if window_size == '1e6' else f'.{window_size}'
    return f'{bucket}/results/misc/ukb.{pop}.for_grm.pruned{cut}.ht'


def get_ukb_grm_plink_path(pop: str, data_iteration: int = 0, window_size: str = '1e6'):
    suffix = f'.{data_iteration}' if data_iteration else ""
    cut = '' if window_size == '1e6' else f'.{window_size}'
    return f'{bucket}/results/misc/ukb.{pop}.for_grm{suffix}.pruned{cut}.plink'


def get_ukb_samples_file_path(pop: str, data_iteration: int = 0):
    suffix = f'.{data_iteration}' if data_iteration else ""
    return f'{bucket}/results/misc/ukb.{pop}{suffix}.samples'