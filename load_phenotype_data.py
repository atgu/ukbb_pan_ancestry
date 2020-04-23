#!/usr/bin/env python3

__author__ = 'konradk'

import argparse
from datetime import date
from ukb_common import *
from ukbb_pan_ancestry.resources import *

pairwise_correlation_ht_path = f'{pheno_folder}/pheno_combo_explore/pairwise_correlations.ht'


def load_dob_ht():
    dob_ht = hl.import_table(pre_phesant_tsv_path, impute=False, min_partitions=100, missing='', key='userId',
                             types={'userId': hl.tint32, 'x54_0_0': hl.tint32})
    year_field, month_field = dob_ht.x34_0_0, dob_ht.x52_0_0
    month_field = hl.cond(hl.len(month_field) == 1, '0' + month_field, month_field)
    dob_ht = dob_ht.select(
        date_of_birth=hl.experimental.strptime(year_field + month_field + '15 00:00:00', '%Y%m%d %H:%M:%S', 'GMT'),
        month_of_birth=month_field,
        year_of_birth=year_field,
        recruitment_center=dob_ht.x54_0_0
    )
    return dob_ht


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
    dob_ht = load_dob_ht()
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
    new_mt = mt.filter_cols(hl.set({'48', '49'}).contains(mt.phenocode) & (mt.coding == 'raw'))
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
    pheno_modifier = hl.if_else(new_mt.phenocode == '48', 'irnt', 'whr')
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


def filter_and_annotate_ukb_data(ht, criteria, type_cast_function = hl.float64):
    fields_to_keep = {x.split('-')[0]: type_cast_function(v) for x, v in ht.row_value.items() if criteria(x, v)}
    ht = ht.select(**fields_to_keep)
    description_ht = load_showcase()
    mt = ht.to_matrix_table_row_major(columns=list(fields_to_keep), entry_field_name='value', col_field_name='pheno')
    return mt.annotate_cols(**description_ht[mt.pheno])


def load_showcase():
    return hl.import_table(pheno_description_path, impute=True, missing='', key='FieldID', types={'FieldID': hl.tstr})


def load_first_occurrence_data(overwrite: bool = False):
    ht = hl.import_table(first_exposure_and_activity_monitor_data_path, delimiter=',', quote='"', missing='', impute=True, key='eid')  #, min_partitions=500)
    pseudo_dates = {'1901-01-01', '2037-07-07'}  # Pseudo date information at http://biobank.ctsu.ox.ac.uk/showcase/coding.cgi?id=819
    dob_ht = load_dob_ht()[ht.key]
    dob = dob_ht.dob
    month = dob_ht.month

    def parse_first_occurrence(x):
        return (hl.case(missing_false=True)
            .when((x.dtype == hl.tint32) | (x.dtype == hl.tfloat64), hl.float64(x))  # Source of the first code ...
            .when(hl.literal(pseudo_dates).contains(hl.str(x)), hl.null(hl.float64))  # Setting past and future dates to missing
            .when(hl.str(x) == '1902-02-02', 0.0)  # Matches DOB
            .when(hl.str(x) == '1903-03-03',  # Within year of birth (taking midpoint between month of birth and EOY)
                  (hl.experimental.strptime('1970-12-31 00:00:00', '%Y-%m-%d %H:%M:%S', 'GMT') -
                   hl.experimental.strptime('1970-' + month + '-15 00:00:00', '%Y-%m-%d %H:%M:%S',
                                            'GMT')) / 2)
            .default(hl.experimental.strptime(hl.str(x) + ' 00:00:00', '%Y-%m-%d %H:%M:%S', 'GMT') - dob
        ))

    mt = filter_and_annotate_ukb_data(ht, lambda k, v: k.startswith('13') and k.endswith('-0.0'), parse_first_occurrence)

    mt = mt.key_cols_by(trait_type='icd_first_occurrence',
                        phenocode=mt.pheno, pheno_sex='both_sexes',
                        coding=NULL_STR_KEY, modifier=NULL_STR_KEY)
    mt.write(get_ukb_pheno_mt_path('icd_first_occurrence'), overwrite)


def load_activity_monitor_data(overwrite: bool = False):
    ht = hl.import_table(first_exposure_and_activity_monitor_data_path, delimiter=',', quote='"', missing='', impute=True, key='eid')  #, min_partitions=500)
    mt = filter_and_annotate_ukb_data(ht, lambda x, v: x.startswith('90') and x.endswith('-0.0') and
                                                       v.dtype in {hl.tint32, hl.tfloat64})
    mt = mt.key_cols_by(trait_type='activity_monitor', pheno=mt.pheno, pheno_sex='both_sexes', modifier=hl.null(hl.tstr), coding=hl.null(hl.tstr))
    mt.write(get_ukb_pheno_mt_path('activity_monitor'), overwrite)


def load_brain_mri_data(overwrite: bool = False):
    ht = hl.import_table(brain_mri_data_path, delimiter=',', quote='"', missing='', impute=True, key='eid')  #, min_partitions=500)
    mt = filter_and_annotate_ukb_data(ht, lambda x, v: v.dtype in {hl.tint32, hl.tfloat64})
    mt = mt.key_cols_by(trait_type='brain_mri', pheno=mt.pheno, pheno_sex='both_sexes', coding=NULL_STR_KEY, modifier=NULL_STR_KEY)
    mt.write(get_ukb_pheno_mt_path('brain_mri'), overwrite)


def load_custom_pheno(traittype_source: str, overwrite: bool = False, sex: str = 'both_sexes'):
    trait_type, source = traittype_source.split('-')
    print(f'Loading {get_custom_pheno_path(traittype_source, extension="txt")}')
    ht = hl.import_table(get_custom_pheno_path(traittype_source, extension='txt'), impute=True)
    inferred_sample_column = list(ht.row)[0]
    ht = ht.key_by(userId=ht[inferred_sample_column]).drop(inferred_sample_column)
    if trait_type == 'categorical':
        ht = ht.annotate(**{x: hl.bool(ht[x]) for x in list(ht.row_value)})

    mt = pheno_ht_to_mt(ht, trait_type, rekey=False).annotate_cols(data_type=trait_type)
    mt = mt.key_cols_by(trait_type=trait_type, phenocode=mt.phesant_pheno, pheno_sex=sex, coding=NULL_STR_KEY,
                        modifier=NULL_STR_KEY).drop('phesant_pheno')
    mt = mt.annotate_cols(category=source)
    return mt.checkpoint(get_custom_pheno_path(traittype_source, extension='ht'), overwrite)


def main(args):
    hl.init(log='/load_pheno.log')
    sexes = ('both_sexes_no_sex_specific', 'females', 'males')
    data_types = ('categorical', 'continuous') #, 'biomarkers')
    curdate = date.today().strftime("%y%m%d")

    if args.load_data:
        load_hesin_data(args.overwrite)  # TODO: create derivative phenotypes from hesin data
        load_first_occurrence_data(args.overwrite)
        load_brain_mri_data(args.overwrite)
        load_activity_monitor_data(args.overwrite)
        load_prescription_data(prescription_tsv_path, prescription_mapping_path).write(get_ukb_pheno_mt_path('prescriptions'), args.overwrite)
        load_icd_data(pre_phesant_tsv_path, icd_codings_ht_path, f'{temp_bucket}/pheno').write(get_ukb_pheno_mt_path('icd'), args.overwrite)
        load_icd_data(pre_phesant_tsv_path, icd9_codings_ht_path, f'{temp_bucket}/pheno_icd9', icd9=True).write(get_ukb_pheno_mt_path('icd9'), args.overwrite)
        icd10 = hl.read_matrix_table(get_ukb_pheno_mt_path('icd')).annotate_cols(icd_version='icd10')
        icd9 = hl.read_matrix_table(get_ukb_pheno_mt_path('icd9')).annotate_cols(icd_version='icd9')
        icd9 = icd9.select_entries('primary_codes', 'secondary_codes', external_codes=hl.null(hl.tbool), cause_of_death_codes=hl.null(hl.tbool), any_codes=icd9.any_codes)
        icd10.union_cols(icd9).write(get_ukb_pheno_mt_path('icd_all'), args.overwrite)

        for sex in sexes:
            ht = hl.import_table(get_phesant_all_phenos_tsv_path(sex), impute=True, min_partitions=100, missing='', key='userId', quote='"', force_bgz=True)
            pheno_ht = ht.checkpoint(get_ukb_pheno_ht_path(sex), overwrite=args.overwrite)
            for data_type in data_types:
                pheno_ht_to_mt(pheno_ht, data_type).write(get_ukb_pheno_mt_path(data_type, sex), args.overwrite)

            ht = hl.import_table(get_ukb_additional_phenos_tsv_path(sex), impute=True, min_partitions=100, missing='NA', key='userId')
            mt = pheno_ht_to_mt(ht, 'continuous', pheno_function_type=hl.str)
            description_ht = hl.import_table(get_ukb_additional_phenos_description_path(), impute=True, quote='"', key=['pheno', 'coding'])
            mt.annotate_cols(**description_ht[mt.col_key]).write(get_ukb_pheno_mt_path('additional', sex), args.overwrite)

        ht = hl.import_table(phesant_biomarker_phenotypes_tsv_path, impute=True, min_partitions=100, missing='', key='userId', force_bgz=True, quote='"')
        ht = ht.checkpoint(get_biomarker_ht_path(), overwrite=args.overwrite)
        mt = pheno_ht_to_mt(ht, 'continuous')
        description_ht = load_showcase()
        mt.annotate_cols(**description_ht[hl.str(mt.pheno)]).write(get_ukb_pheno_mt_path('biomarkers'), args.overwrite)

    if args.combine_data:
        for data_type in data_types:
            mt = combine_datasets({sex: get_ukb_pheno_mt_path(data_type, sex) for sex in sexes},
                                  {sex: get_ukb_phesant_summary_tsv_path(sex) for sex in sexes},
                                  pheno_description_path, coding_ht_path, data_type)
            mt.write(get_ukb_pheno_mt_path(data_type, 'full'), args.overwrite)

        mt = combine_datasets({sex: get_ukb_pheno_mt_path('additional', sex) for sex in sexes})
        mt.write(get_ukb_pheno_mt_path('additional', 'full'), args.overwrite)

        data_types = ('categorical', 'continuous', 'biomarkers', 'icd_all', 'prescriptions', 'phecode', 'additional',
                      'brain_mri', 'activity_monitor', 'icd_first_occurrence')
        location = {'categorical': 'full', 'continuous': 'full', 'additional': 'full'}
        pheno_file_dict = {data_type: hl.read_matrix_table(get_ukb_pheno_mt_path(data_type, location.get(data_type, 'both_sexes_no_sex_specific')))
                           for data_type in data_types}
        cov_ht = get_covariates(hl.int32).persist()
        mt = combine_pheno_files_multi_sex(pheno_file_dict, cov_ht)
        mt = add_white_noise_pheno(mt)
        mt = mt.checkpoint(f'{temp_bucket}/pheno_{curdate}.mt', overwrite=args.overwrite, _read_if_exists=not args.overwrite)
        mt = add_whr(mt)
        mt = mt.checkpoint(get_ukb_pheno_mt_path(), args.overwrite, _read_if_exists=not args.overwrite)
    summarize_data(args.overwrite)


    if args.add_dataset:
        mt = load_custom_pheno(args.add_dataset, args.overwrite)
        cov_ht = get_covariates(hl.int32).persist()
        mt = combine_pheno_files_multi_sex({'custom': mt}, cov_ht)
        original_mt = hl.read_matrix_table(get_ukb_pheno_mt_path())
        original_mt = original_mt.checkpoint(get_ukb_pheno_mt_path(f'full_before_{curdate}', sex='full'), args.overwrite)
        original_mt.cols().export(f'{pheno_folder}/all_pheno_summary_before_{curdate}.txt.bgz')
        mt = original_mt.union_cols(mt).write(get_ukb_pheno_mt_path(), args.overwrite)
        summarize_data(args.overwrite)


    if args.summarize_data:
        summarize_data(args.overwrite)
        ht = hl.read_table(get_phenotype_summary_path('full'))

        ht = ht.filter(
            ~hl.set({'raw', 'icd9'}).contains(ht.coding) &
            ~hl.set({'22601', '22617', '20024', '41230', '41210'}).contains(ht.pheno) &
            (ht.n_cases_both_sexes >= MIN_CASES_ALL) &
            (ht.n_cases_by_pop >= hl.cond(ht.pop == 'EUR', MIN_CASES_EUR, MIN_CASES))
        )
        ht = ht.group_by('pop').aggregate(n_phenos=hl.agg.count(),
                                          n_samples=hl.agg.max(ht.n_cases_by_pop))
        ht.show()
        print(f'Total of {ht.aggregate(hl.agg.sum(ht.n_phenos))} phenos')

    if args.pairwise_correlations:
        mt = hl.read_matrix_table(get_ukb_pheno_mt_path())
        make_correlation_ht(mt).write(pairwise_correlation_ht_path, args.overwrite)
        hl.read_table(pairwise_correlation_ht_path).flatten().export(pairwise_correlation_ht_path.replace('.ht', '.txt.bgz'))


def summarize_data(overwrite):
    mt = hl.read_matrix_table(get_ukb_pheno_mt_path())
    ht = mt.group_rows_by('pop').aggregate(
        stats=hl.agg.stats(mt.both_sexes),
        n_cases_by_pop=hl.cond(hl.set({'continuous', 'biomarkers'}).contains(mt.data_type),
                               hl.agg.count_where(hl.is_defined(mt.both_sexes)),
                               hl.int64(hl.agg.sum(mt.both_sexes)))
    ).entries()
    ht = ht.key_by('pop', 'pheno', 'coding', trait_type=ht.data_type)
    ht = ht.checkpoint(get_phenotype_summary_path('full'), overwrite=overwrite, _read_if_exists=not overwrite)
    ht.flatten().export(get_phenotype_summary_path('full', 'tsv'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--load_data', help='Load data', action='store_true')
    parser.add_argument('--combine_data', help='Load data', action='store_true')
    parser.add_argument('--add_dataset', help='Load data')
    parser.add_argument('--pairwise_correlations', help='Load data', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
