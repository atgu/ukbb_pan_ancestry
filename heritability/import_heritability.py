__author__ = 'Rahul Gupta'

import hail as hl
from ukbb_common.resources.generic import PHENO_KEY_FIELDS


def get_h2_flat_file():
    return 'gs://ukb-diverse-pops-public-free/h2/h2_estimates_all_flat_220317.tsv'
    #return 'gs://ukb-diverse-pops/rg-h2-tables/h2_estimates_all_flat_211101.tsv'


def get_h2_ht():
    return 'gs://ukb-diverse-pops-public-free/h2/h2_estimates_all.ht'
    #return 'gs://ukb-diverse-pops/rg-h2-tables/ht/h2_estimates_all.ht'


def qc_to_flags(qc_struct):
    map_flags_opposite = {'GWAS_run': 'GWAS_not_run',
                          'ancestry_reasonable_n': 'ancestry_n_too_low',
                          'defined_h2': 'h2_not_defined',
                          'significant_z': 'h2_z_insignificant',
                          'in_bounds_h2': 'out_of_bounds_h2',
                          'normal_lambda': 'out_of_bounds_lambda',
                          'normal_ratio': 'fail_ratio',
                          'EUR_plus_1': 'not_EUR_plus_1'}
    qc_struct = qc_struct.rename(map_flags_opposite)
    qc_flags = [x for x in qc_struct.keys() if x not in ["pass_all"]]
    first_failed_qc = hl.zip(qc_flags, [~qc_struct[x] for x in qc_flags]).filter(lambda x: x[1]).first()[0]
    return hl.if_else(qc_struct.pass_all, 'PASS', first_failed_qc)


def import_h2_flat_file(save_to_ht, overwrite):
    ht = hl.import_table(get_h2_flat_file(), 
                         delimiter='\t', 
                         impute=True, 
                         key=PHENO_KEY_FIELDS)
    ht = ht.rename({'ancestry':'pop'})
    ht = ht.drop('phenotype_id')

    # solution from Zulip to munge flat file columns into nested structs
    def recur(dict_ref, split_name):
        if (len(split_name) == 1):
            dict_ref[split_name[0]] = row[name]
            return
        existing = dict_ref.get(split_name[0])
        if existing is not None:
            assert isinstance(existing, dict), existing
            recur(existing, split_name[1:])
        else:
            existing = {}
            dict_ref[split_name[0]] = existing
            recur(existing, split_name[1:])


    d = {}
    row = ht.row
    for name in row:
        if name not in ht.key:
            recur(d, name.split('.'))
    
    
    def dict_to_struct(d):
        fields = {}
        for k, v in d.items():
            if isinstance(v, dict):
                v = dict_to_struct(v)
            fields[k] = v
        return hl.struct(**fields)
    
    
    ht = ht.select(**dict_to_struct(d))
    ht_collect = ht.collect_by_key().rename({'values':'heritability'})
    ht_collect_sorted = ht_collect.annotate(heritability = hl.sorted(ht_collect.heritability, key=lambda x: x.pop))
    ht_collect_sorted = ht_collect_sorted.repartition(250)

    if save_to_ht:
        ht_collect_sorted.write(get_h2_ht(), overwrite=overwrite)

    return ht_collect_sorted
