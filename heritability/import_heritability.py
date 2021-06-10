__author__ = 'Rahul Gupta'

import hail as hl


def get_h2_flat_file():
    return 'gs://ukb-diverse-pops/rg-h2-tables/h2_estimates_all_flat.tsv'


def get_h2_ht():
    return 'gs://ukb-diverse-pops/rg-h2-tables/ht/h2_estimates_all.ht'


def import_h2_flat_file(save_to_ht, overwrite):
    ht = hl.import_table(get_h2_flat_file(), 
                         delimiter='\t', 
                         impute=True, 
                         key=['trait_type','phenocode','pheno_sex','coding','modifier'])
    ht = ht.rename({'ancestry':'pop'})
    ht = ht.drop('phenotype_id')
    # solution from Zulip to munge flat file columns into nested structs
    d = {}
    row = ht.row
    for name in row:
        if name not in ht.key:
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
