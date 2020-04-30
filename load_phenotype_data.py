#!/usr/bin/env python3

__author__ = 'konradk'

import argparse
from datetime import date
from ukbb_pan_ancestry import *


def main(args):
    hl.init(log='/load_pheno.log')
    sexes = ('both_sexes_no_sex_specific', 'females', 'males')
    data_types = ('categorical', 'continuous') #, 'biomarkers')
    curdate = date.today().strftime("%y%m%d")

    # load_brain_mri_data(brain_mri_data_path).write(get_ukb_pheno_mt_path('brain_mri'), overwrite)
    if args.load_data:
        load_first_occurrence_data(first_exposure_and_activity_monitor_data_path,
                                   pre_phesant_tsv_path).write(get_ukb_pheno_mt_path('icd_first_occurrence'), args.overwrite)
        load_activity_monitor_data(first_exposure_and_activity_monitor_data_path).write(get_ukb_pheno_mt_path('activity_monitor'), args.overwrite)
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
                      'activity_monitor', 'icd_first_occurrence')  # 'brain_mri')
        location = {'categorical': 'full', 'continuous': 'full', 'additional': 'full'}
        pheno_file_dict = {data_type: hl.read_matrix_table(get_ukb_pheno_mt_path(data_type, location.get(data_type, 'both_sexes_no_sex_specific')))
                           for data_type in data_types}
        cov_ht = get_covariates(hl.int32).persist()
        mt = combine_pheno_files_multi_sex(pheno_file_dict, cov_ht)
        mt = add_white_noise_pheno(mt)
        mt = mt.checkpoint(f'{temp_bucket}/pheno_{curdate}.mt', overwrite=args.overwrite, _read_if_exists=not args.overwrite)
        mt = add_whr(mt)
        mt.write(get_ukb_pheno_mt_path(), args.overwrite)
        summarize_data(args.overwrite)

    if args.add_dataset:
        if args.add_dataset == 'covid':
            mt = load_covid_data(pre_phesant_tsv_path).checkpoint(get_ukb_pheno_mt_path('covid'), args.overwrite)
        else:
            mt = load_custom_pheno(args.add_dataset).checkpoint(get_custom_pheno_path(args.add_dataset, extension='ht'), args.overwrite)
        cov_ht = get_covariates(hl.int32).persist()
        mt = combine_pheno_files_multi_sex({'custom': mt}, cov_ht)
        original_mt = hl.read_matrix_table(get_ukb_pheno_mt_path())
        original_mt = original_mt.checkpoint(get_ukb_pheno_mt_path(f'full_before_{curdate}', sex='full'), args.overwrite)
        original_mt.cols().export(f'{pheno_folder}/all_pheno_summary_before_{curdate}.txt.bgz')
        original_mt.union_cols(mt, row_join_type='outer').write(get_ukb_pheno_mt_path(), args.overwrite)
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
        make_pairwise_ht(mt, mt.both_sexes, correlation=True).write(pairwise_correlation_ht_path, args.overwrite)
        hl.read_table(pairwise_correlation_ht_path).flatten().export(pairwise_correlation_ht_path.replace('.ht', '.txt.bgz'))

        mt = hl.read_matrix_table(get_ukb_pheno_mt_path())
        bool_types = {'icd10', 'prescriptions', 'categorical'}
        mt = mt.select_entries(pheno_indicator=hl.case(missing_false=True)
                               .when((mt.trait_type == 'icd_first_occurrence') & mt.description.startswith('Date'), hl.is_defined(mt.both_sexes))
                               .when(hl.literal(bool_types).contains(mt.trait_type), hl.bool(mt.both_sexes))
                               .default(hl.null(hl.tbool)))
        mt = mt.filter_cols(hl.agg.all(hl.is_defined(mt.pheno_indicator)))
        new_mt = combine_phenotypes_in_mt(mt, {'combined_parkinsons': ["dopamine pro-drug _ decarboxylase inhibitor|Parkinson's", "131022"]})

        mt = mt.union_cols(new_mt)
        ht = make_pairwise_ht(mt, mt.pheno_indicator)
        ht.write(pairwise_cooccurrence_ht_path, True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--load_data', help='Load data', action='store_true')
    parser.add_argument('--combine_data', help='Load data', action='store_true')
    parser.add_argument('--summarize_data', help='Load data', action='store_true')
    parser.add_argument('--add_dataset', help='Load data')
    parser.add_argument('--pairwise_correlations', help='Load data', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        from slack_token_pkg.slack_creds import slack_token
        with slack.slack_notifications(slack_token, args.slack_channel):
            main(args)
