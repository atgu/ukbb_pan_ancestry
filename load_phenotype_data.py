#!/usr/bin/env python3

__author__ = 'konradk'

from gnomad_hail import *
from ukb_phenotypes import *
from ukbb_pan_ancestry.resources import *

temp_bucket = 'gs://ukbb-diverse-temp-30day'


def main(args):
    hl.init(log='/load_pheno.log')
    # sexes = ('both_sexes_no_sex_specific', 'females', 'males')
    data_types = ('categorical', 'continuous') #, 'biomarkers')

    if args.load_data:
        load_icd_data(pre_phesant_tsv_path, icd_codings_path, f'{temp_bucket}/pheno').write(get_ukb_pheno_mt_path('icd'), args.overwrite)
        load_icd_data(pre_phesant_tsv_path, icd9_codings_path, f'{temp_bucket}/pheno_icd9', icd9=True).write(get_ukb_pheno_mt_path('icd9'), args.overwrite)
        icd10 = hl.read_matrix_table(get_ukb_pheno_mt_path('icd')).annotate_cols(icd_version='icd10')
        icd9 = hl.read_matrix_table(get_ukb_pheno_mt_path('icd9')).annotate_cols(icd_version='icd9')
        icd9 = icd9.select_entries('primary_codes', 'secondary_codes', external_codes=hl.null(hl.tbool), cause_of_death_codes=hl.null(hl.tbool), any_codes=icd9.any_codes)
        icd10.union_cols(icd9).write(get_ukb_pheno_mt_path('icd_all'), args.overwrite)
        # read_covariate_data(get_pre_phesant_data_path()).write(get_ukb_covariates_ht_path(), args.overwrite)

        # for sex in sexes:
        # ht = hl.import_table(phesant_all_phenos_tsv_paths.format(1), impute=True, min_partitions=100, missing='', key='userId', quote='"', force_bgz=True)
        # for i in range(2, 4):
        #     pheno_ht = hl.import_table(phesant_all_phenos_tsv_paths.format(i), impute=True, min_partitions=100, missing='', key='userId', quote='"', force_bgz=True)
        #     ht = ht.annotate(**pheno_ht[ht.key])
        # ht.write(get_ukb_pheno_ht_path(), overwrite=args.overwrite)

        pheno_ht = hl.read_table(get_ukb_pheno_ht_path())
        for data_type in data_types:
            pheno_ht_to_mt(pheno_ht, data_type).write(get_ukb_pheno_mt_path(data_type), args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--load_data', help='Load data', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
