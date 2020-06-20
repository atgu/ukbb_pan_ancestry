from ukb_common import *
from ukbb_pan_ancestry.resources import *


def add_white_noise_pheno(mt):
    new_mt = mt.add_col_index()
    new_mt = (new_mt.filter_cols(new_mt.col_idx == 0).drop('col_idx')
              .key_cols_by(trait_type='continuous', phenocode='random', pheno_sex='both_sexes', coding='', modifier='random'))
    new_mt = new_mt.select_entries(both_sexes=hl.rand_norm(seed=42),
                                   females=hl.or_missing(new_mt.sex == 0, hl.rand_norm(seed=43)),
                                   males=hl.or_missing(new_mt.sex == 1, hl.rand_norm(mean=1, seed=44)))
    new_mt = new_mt.select_cols(n_cases_both_sexes=hl.agg.count(),
                                n_cases_females=hl.agg.count_where(new_mt.sex == 0),
                                n_cases_males=hl.agg.count_where(new_mt.sex == 1),
                                description='hl.rand_norm(seed=42)', description_more=NULL_STR,
                                coding_description=NULL_STR, category=NULL_STR)

    new_mt2 = mt.add_col_index()
    pop_dict = new_mt2.aggregate_rows(
        hl.dict(hl.zip_with_index(hl.array(hl.agg.collect_as_set(new_mt2.pop)), index_first=False)),
        _localize=False)
    new_mt2 = (new_mt2.filter_cols(new_mt2.col_idx == 0).drop('col_idx')
              .key_cols_by(trait_type='continuous', phenocode='random', pheno_sex='both_sexes', coding='', modifier='random_strat'))
    new_mt2 = new_mt2.select_entries(both_sexes=hl.rand_norm(mean=pop_dict[new_mt2.pop], seed=42),
                                     females=hl.or_missing(new_mt2.sex == 0, hl.rand_norm(mean=pop_dict[new_mt2.pop], seed=43)),
                                     males=hl.or_missing(new_mt2.sex == 1, hl.rand_norm(mean=pop_dict[new_mt2.pop] + 1, seed=44)))
    new_mt2 = new_mt2.select_cols(n_cases_both_sexes=hl.agg.count(),
                                  n_cases_females=hl.agg.count_where(new_mt2.sex == 0),
                                  n_cases_males=hl.agg.count_where(new_mt2.sex == 1),
                                  description='hl.rand_norm(mean=pop_dict[mt.pop], seed=42)', description_more=NULL_STR,
                                  coding_description=NULL_STR, category=NULL_STR)

    return mt.union_cols(new_mt).union_cols(new_mt2)


def irnt(he: hl.expr.Expression, output_loc: str = 'irnt'):
    ht = he._indices.source
    n_rows = ht.count()
    ht = ht.order_by(he).add_index()
    return ht.annotate(**{output_loc: hl.qnorm((ht.idx + 0.5) / n_rows)})


def add_whr(mt):
    new_mt = mt.filter_cols(hl.set({'48', '49'}).contains(mt.phenocode) & (mt.modifier == 'raw'))
    new_ht = new_mt.annotate_rows(whr=hl.agg.sum(new_mt.both_sexes * hl.int(new_mt.phenocode == '48')) /
                                      hl.agg.sum(new_mt.both_sexes * hl.int(new_mt.phenocode == '49'))
                                  ).rows()
    new_ht = irnt(new_ht.whr).key_by('userId').select('whr', 'irnt')
    new_mt = new_mt.annotate_rows(**new_ht[new_mt.row_key])
    # Hijacking the 48 and 49 phenos to be irnt and raw respectively
    pheno_value = hl.if_else(new_mt.phenocode == '48', new_mt.irnt, new_mt.whr)
    new_mt = new_mt.annotate_entries(both_sexes=pheno_value,
                                     females=hl.or_missing(new_mt.sex == 0, pheno_value),
                                     males=hl.or_missing(new_mt.sex == 1, pheno_value))
    pheno_modifier = hl.if_else(new_mt.phenocode == '48', 'irnt', 'raw')
    new_mt = (new_mt
              .key_cols_by(trait_type='continuous',
                           phenocode='whr', pheno_sex='both_sexes',
                           coding=NULL_STR_KEY,
                           modifier=pheno_modifier))
    new_mt = new_mt.select_cols(n_cases_both_sexes=hl.agg.count_where(hl.is_defined(new_mt.both_sexes)),
                                n_cases_females=hl.agg.count_where(hl.is_defined(new_mt.females)),
                                n_cases_males=hl.agg.count_where(hl.is_defined(new_mt.males)),
                                description='pheno 48 / pheno 49', description_more=NULL_STR,
                                coding_description=NULL_STR, category=NULL_STR)

    return mt.union_cols(new_mt)


def load_brain_mri_data(brain_mri_data_path):
    ht = hl.import_table(brain_mri_data_path, delimiter=',', quote='"', missing='', impute=True, key='eid')  #, min_partitions=500)
    mt = filter_and_annotate_ukb_data(ht, lambda x, v: v.dtype in {hl.tint32, hl.tfloat64})
    mt = mt.key_cols_by(trait_type='brain_mri', phenocode=mt.phenocode, pheno_sex='both_sexes', coding=NULL_STR_KEY, modifier=NULL_STR_KEY)
    return mt



def load_hesin_data(overwrite: bool = False):
    overall_ht = hl.import_table(get_hesin_raw_data_path(), impute=True, min_partitions=40, missing='', key=('eid', 'ins_index'))
    overall_ht = overall_ht.filter(overall_ht.dsource == 'HES')

    # Useful fields:
    # pctcode	Current PCT responsible for patient
    # gpprpct	PCT where patients GP was registered
    # More coding data in:
    # coding_ht = hl.read_table(coding_ht_path)

    diag_ht = hl.import_table(get_hesin_raw_data_path('diag'), impute=True, min_partitions=40, missing='', key=('eid', 'ins_index'))
    diag_ht = diag_ht.annotate(**overall_ht[diag_ht.key])
    date_fields = ('epistart', 'epiend', 'elecdate', 'admidate', 'disdate')
    # diag_ht = diag_ht.annotate(admidate=hl.experimental.strptime(diag_ht.admidate, '%Y%m%d'))
    diag_ht = diag_ht.annotate(**{date: hl.experimental.strptime(hl.str(diag_ht[date]) + ' 00:00:00', '%Y%m%d %H:%M:%S', 'GMT') for date in date_fields})
    # assert diag_ht.aggregate(hl.agg.all(hl.is_defined(diag_ht.admidate)))
    diag_ht = diag_ht.key_by('eid', 'diag_icd10').collect_by_key()
    # Only first entry by admission date
    diag_ht = diag_ht.annotate(values=hl.sorted(diag_ht.values, key=lambda x: x.admidate)[0])
    diag_ht = diag_ht.transmute(**diag_ht.values)
    diag_mt = diag_ht.to_matrix_table(row_key=['eid'], col_key=['diag_icd10'])
    dob_ht = load_dob_ht(pre_phesant_tsv_path)
    diag_mt = diag_mt.annotate_rows(dob=dob_ht[diag_mt.row_key].dob)
    diag_mt.write(get_hesin_mt_path('diag'), overwrite)

    oper_ht = hl.import_table(get_hesin_raw_data_path('oper'), impute=True, missing='', key=('eid', 'ins_index'))
    oper_ht = oper_ht.annotate(**overall_ht[oper_ht.key])
    oper_ht = oper_ht.filter(hl.is_defined(oper_ht.admidate))
    oper_ht = oper_ht.key_by('eid', 'oper4').collect_by_key('operations')
    oper_mt = oper_ht.to_matrix_table(row_key=['eid'], col_key=['oper4'])
    oper_mt.write(get_hesin_mt_path('oper'), overwrite)
    # Useful fields:
    # posopdur	Duration of post-operative stay

    maternity_ht = hl.import_table(get_hesin_raw_data_path('maternity'), impute=True, missing='', key=('eid', 'ins_index'))
    delivery_ht = hl.import_table(get_hesin_raw_data_path('delivery'), impute=True, missing='', key=('eid', 'ins_index'))
    delivery_ht = delivery_ht.annotate(**maternity_ht[delivery_ht.key])
    delivery_ht = delivery_ht.annotate(**overall_ht[delivery_ht.key])
    delivery_ht.write(get_hesin_delivery_ht_path(), overwrite)


def load_custom_pheno(traittype_source: str, sex: str = 'both_sexes', extension: str = 'txt'):
    trait_type, source = traittype_source.split('-')
    print(f'Loading {get_custom_pheno_path(traittype_source, extension=extension)}')
    ht = hl.import_table(get_custom_pheno_path(traittype_source, extension=extension), impute=True)
    inferred_sample_column = list(ht.row)[0]
    ht = ht.key_by(userId=ht[inferred_sample_column])
    if inferred_sample_column != 'userId':
        ht = ht.drop(inferred_sample_column)
    if trait_type == 'categorical':
        ht = ht.annotate(**{x: hl.bool(ht[x]) for x in list(ht.row_value)})

    mt = pheno_ht_to_mt(ht, trait_type, rekey=False).annotate_cols(data_type=trait_type)
    mt = mt.key_cols_by(trait_type=trait_type, phenocode=mt.phesant_pheno, pheno_sex=sex, coding=NULL_STR_KEY,
                        modifier=NULL_STR_KEY).drop('phesant_pheno')
    mt = mt.annotate_cols(category=source)
    return mt


def summarize_data(overwrite):
    mt = hl.read_matrix_table(get_ukb_pheno_mt_path())
    ht = mt.group_rows_by('pop').aggregate(
        stats=hl.agg.stats(mt.both_sexes),
        n_cases_by_pop=hl.cond(hl.set({'continuous', 'biomarkers'}).contains(mt.trait_type),
                               hl.agg.count_where(hl.is_defined(mt.both_sexes)),
                               hl.int64(hl.agg.sum(mt.both_sexes)))
    ).entries()
    ht = ht.key_by('pop', *PHENO_KEY_FIELDS)
    ht = ht.checkpoint(get_phenotype_summary_path('full'), overwrite=overwrite, _read_if_exists=not overwrite)
    ht.flatten().export(get_phenotype_summary_path('full', 'tsv'))


def combine_phenotypes_in_mt(mt, pheno_combo_dict):
    new_mt = combine_phenotypes_with_name(mt.key_cols_by(), mt.phenocode, mt.pheno_indicator,
                                          pheno_combo_dict, 'new_pheno', 'pheno_indicator')
    new_mt = new_mt.key_cols_by(
        **{x: new_mt.new_pheno if x == 'phenocode' else '' for x in PHENO_KEY_FIELDS})
    new_mt = new_mt.annotate_cols(
        **compute_cases_binary(new_mt.pheno_indicator, new_mt.sex),
        **{x: hl.null(hl.tstr) for x in PHENO_COLUMN_FIELDS if 'n_cases' not in x})
    return new_mt
