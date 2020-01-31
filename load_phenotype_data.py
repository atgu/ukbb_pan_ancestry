#!/usr/bin/env python3

__author__ = 'konradk'

from gnomad_hail import *
from ukb_common import *
from ukbb_pan_ancestry.resources import *

temp_bucket = 'gs://ukbb-diverse-temp-30day'
pairwise_correlation_ht_path = f'{pheno_folder}/pheno_combo_explore/pairwise_correlations.ht'


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
    dob_ht = hl.import_table(pre_phesant_tsv_path, impute=False, min_partitions=100, missing='', key='userId', types={'userId': hl.tint32})
    year_field, month_field = dob_ht.x34_0_0, dob_ht.x52_0_0
    month_field = hl.cond(hl.len(month_field) == 1, '0' + month_field, month_field)
    dob_ht = dob_ht.select(dob=hl.experimental.strptime(year_field + month_field + '15 00:00:00', '%Y%m%d %H:%M:%S', 'GMT'))
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
              .key_cols_by(pheno='random', coding='random')
              .annotate_cols(data_type='continuous', meaning='hl.rand_norm(seed=42)', path=''))
    new_mt = new_mt.annotate_entries(both_sexes=hl.rand_norm(seed=42),
                                     females=hl.or_missing(new_mt.sex == 0, hl.rand_norm(seed=43)),
                                     males=hl.or_missing(new_mt.sex == 1, hl.rand_norm(mean=1, seed=44)))
    new_mt = new_mt.annotate_cols(n_cases_both_sexes=hl.agg.count(),
                                  n_cases_females=hl.agg.count_where(new_mt.sex == 0),
                                  n_cases_males=hl.agg.count_where(new_mt.sex == 1))

    new_mt2 = mt.add_col_index()
    pop_dict = new_mt2.aggregate_rows(
        hl.dict(hl.zip_with_index(hl.array(hl.agg.collect_as_set(new_mt2.pop)), index_first=False)),
        _localize=False)
    new_mt2 = (new_mt2.filter_cols(new_mt2.col_idx == 0).drop('col_idx')
              .key_cols_by(pheno='random', coding='random_strat')
              .annotate_cols(data_type='continuous', meaning='hl.rand_norm(seed=42)', path=''))
    new_mt2 = new_mt2.annotate_entries(both_sexes=hl.rand_norm(mean=pop_dict[new_mt2.pop], seed=42),
                                     females=hl.or_missing(new_mt2.sex == 0, hl.rand_norm(mean=pop_dict[new_mt2.pop], seed=43)),
                                     males=hl.or_missing(new_mt2.sex == 1, hl.rand_norm(mean=pop_dict[new_mt2.pop] + 1, seed=44)))
    new_mt2 = new_mt2.annotate_cols(n_cases_both_sexes=hl.agg.count(),
                                  n_cases_females=hl.agg.count_where(new_mt2.sex == 0),
                                  n_cases_males=hl.agg.count_where(new_mt2.sex == 1))

    return mt.union_cols(new_mt).union_cols(new_mt2)


def irnt(he: hl.expr.Expression, output_loc: str = 'irnt'):
    ht = he._indices.source
    n_rows = ht.count()
    ht = ht.order_by(he).add_index()
    return ht.annotate(**{output_loc: hl.qnorm((ht.idx + 0.5) / n_rows)})


def add_whr(mt):
    new_mt = mt.filter_cols(hl.set({'48', '49'}).contains(mt.pheno) & (mt.coding == 'raw'))
    new_ht = new_mt.annotate_rows(whr=hl.agg.sum(new_mt.both_sexes * hl.int(new_mt.pheno == '48')) /
                                      hl.agg.sum(new_mt.both_sexes * hl.int(new_mt.pheno == '49'))
                                  ).rows()
    new_ht = irnt(new_ht.whr).key_by('userId')
    new_mt = new_mt.annotate_rows(irnt=new_ht[new_mt.row_key].irnt)
    new_mt = new_mt.filter_cols(new_mt.pheno == '48')
    new_mt = new_mt.annotate_entries(both_sexes=new_mt.irnt,
                                     females=hl.or_missing(new_mt.sex == 0, new_mt.irnt),
                                     males=hl.or_missing(new_mt.sex == 1, new_mt.irnt))
    new_mt = (new_mt
              .key_cols_by(pheno='whr', coding='whr')
              .annotate_cols(data_type='continuous', meaning='pheno 48 / pheno 49', path=''))
    new_mt = new_mt.annotate_cols(n_cases_both_sexes=hl.agg.count_where(hl.is_defined(new_mt.both_sexes)),
                                  n_cases_females=hl.agg.count_where(hl.is_defined(new_mt.females)),
                                  n_cases_males=hl.agg.count_where(hl.is_defined(new_mt.males)))

    return mt.union_cols(new_mt)


def main(args):
    hl.init(log='/load_pheno.log')
    sexes = ('both_sexes_no_sex_specific', 'females', 'males')
    data_types = ('categorical', 'continuous') #, 'biomarkers')

    if args.load_data:
        load_hesin_data(args.overwrite)
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
            pheno_ht_to_mt(ht, 'continuous', pheno_function_type=hl.str).write(get_ukb_pheno_mt_path('additional', sex), args.overwrite)

        ht = hl.import_table(phesant_biomarker_phenotypes_tsv_path, impute=True, min_partitions=100, missing='', key='userId', force_bgz=True, quote='"')
        ht = ht.checkpoint(get_biomarker_ht_path(), overwrite=args.overwrite)
        pheno_ht_to_mt(ht, 'continuous').write(get_ukb_pheno_mt_path('biomarkers'), args.overwrite)

    if args.combine_data:
        for data_type in data_types:
            mt = combine_datasets({sex: get_ukb_pheno_mt_path(data_type, sex) for sex in sexes},
                                  {sex: get_ukb_phesant_summary_tsv_path(sex) for sex in sexes},
                                  pheno_description_path, coding_ht_path, data_type)
            mt.write(get_ukb_pheno_mt_path(data_type, 'full'), args.overwrite)

        mt = combine_datasets({sex: get_ukb_pheno_mt_path('additional', sex) for sex in sexes})
        mt.write(get_ukb_pheno_mt_path('additional', 'full'), args.overwrite)

        data_types = ('categorical', 'continuous', 'biomarkers', 'icd_all', 'prescriptions', 'phecode', 'additional')
        location = {'categorical': 'full', 'continuous': 'full', 'additional': 'full'}
        pheno_file_dict = {data_type: hl.read_matrix_table(get_ukb_pheno_mt_path(data_type, location.get(data_type, 'both_sexes_no_sex_specific')))
                           for data_type in data_types}
        cov_ht = get_covariates(hl.int32).persist()
        mt = combine_pheno_files_multi_sex(pheno_file_dict, cov_ht)
        mt = add_white_noise_pheno(mt)
        mt = add_whr(mt)
        mt = mt.checkpoint(get_ukb_pheno_mt_path(), args.overwrite, _read_if_exists=not args.overwrite)
        mt.cols().export(f'{pheno_folder}/all_pheno_summary.txt.bgz')

        ht = mt.group_rows_by('pop').aggregate(
            stats=hl.agg.stats(mt.both_sexes),
            n_cases_by_pop=hl.cond(hl.set({'continuous', 'biomarkers'}).contains(mt.data_type),
                                   hl.agg.count_where(hl.is_defined(mt.both_sexes)),
                                   hl.int64(hl.agg.sum(mt.both_sexes)))
        ).entries()
        ht = ht.key_by('pop', 'pheno', 'coding', trait_type=ht.data_type)
        ht = ht.checkpoint(get_phenotype_summary_path('full'), overwrite=args.overwrite, _read_if_exists=not args.overwrite)
        ht.flatten().export(get_phenotype_summary_path('full', 'tsv'))

        ht = ht.filter(
            ~hl.set({'raw', 'icd9'}).contains(ht.coding) &
            ~hl.set({'22601', '22617', '20024', '41230', '41210'}).contains(ht.pheno) &
            (ht.n_cases_both_sexes >= MIN_CASES_ALL) &
            (ht.n_cases_by_pop >= hl.cond(ht.pop == 'EUR', MIN_CASES_EUR, MIN_CASES))
        )
        ht = ht.group_by('pop').aggregate(n_phenos=hl.agg.count())
        ht.show()
        print(f'Total of {ht.aggregate(hl.agg.sum(ht.n_phenos))} phenos')

    if args.pairwise_correlations:
        mt = hl.read_matrix_table(get_ukb_pheno_mt_path())
        make_correlation_ht(mt).write(pairwise_correlation_ht_path, args.overwrite)
        hl.read_table(pairwise_correlation_ht_path).flatten().export(pairwise_correlation_ht_path.replace('.ht', '.txt.bgz'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--load_data', help='Load data', action='store_true')
    parser.add_argument('--combine_data', help='Load data', action='store_true')
    parser.add_argument('--pairwise_correlations', help='Load data', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
