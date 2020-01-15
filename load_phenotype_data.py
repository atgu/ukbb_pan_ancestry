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
        # read_covariate_data(get_pre_phesant_data_path()).write(get_ukb_covariates_ht_path(), args.overwrite)

        for sex in sexes:
            ht = hl.import_table(get_phesant_all_phenos_tsv_path(sex), impute=True, min_partitions=100, missing='', key='userId', quote='"', force_bgz=True)
            pheno_ht = ht.checkpoint(get_ukb_pheno_ht_path(sex), overwrite=args.overwrite)
            for data_type in data_types:
                pheno_ht_to_mt(pheno_ht, data_type).write(get_ukb_pheno_mt_path(data_type, sex), args.overwrite)

        ht = hl.import_table(phesant_biomarker_phenotypes_tsv_path, impute=True, min_partitions=100, missing='', key='userId', force_bgz=True, quote='"')
        ht = ht.checkpoint(get_biomarker_ht_path(), overwrite=args.overwrite)
        pheno_ht_to_mt(ht, 'continuous').write(get_ukb_pheno_mt_path('biomarkers'), args.overwrite)

    if args.combine_data:
        for data_type in data_types:
            mt = combine_datasets({sex: get_ukb_pheno_mt_path(data_type, sex) for sex in sexes},
                                  {sex: get_ukb_phesant_summary_tsv_path(sex) for sex in sexes},
                                  pheno_description_path, coding_ht_path, data_type)
            mt.write(get_ukb_pheno_mt_path(data_type, 'full'), args.overwrite)

        data_types = ('categorical', 'continuous', 'biomarkers', 'icd_all', 'prescriptions', 'phecode')
        location = {'categorical': 'full', 'continuous': 'full'}
        pheno_file_dict = {data_type: hl.read_matrix_table(get_ukb_pheno_mt_path(data_type, location.get(data_type, 'both_sexes_no_sex_specific')))
                           for data_type in data_types}
        cov_ht = hl.import_table(get_ukb_meta_pop_tsv_path(), key='s', impute=True)
        cov_ht = cov_ht.annotate(sex=pheno_file_dict['biomarkers'].rows()[cov_ht.key].sex)
        combine_pheno_files_multi_sex(pheno_file_dict, cov_ht).write(get_ukb_pheno_mt_path(), args.overwrite)
        hl.read_matrix_table(get_ukb_pheno_mt_path()).cols().export(f'{pheno_folder}/all_pheno_summary.txt.bgz')

        mt = hl.read_matrix_table(get_ukb_pheno_mt_path())
        make_correlation_ht(mt).write(pairwise_correlation_ht_path, args.overwrite)
        hl.read_table(pairwise_correlation_ht_path).flatten().export(pairwise_correlation_ht_path.replace('.ht', '.txt.bgz'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--load_data', help='Load data', action='store_true')
    parser.add_argument('--combine_data', help='Load data', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
