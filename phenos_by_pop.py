from gnomad_hail import *
from ukb_common import *
from ukbb_pan_ancestry import *

temp_bucket = 'gs://ukbb-diverse-temp-30day'
total_samples = 487409


def main():
    hl.init(default_reference='GRCh37', log='/pheno_agg.log')

    if args.pheno_summary:
        mt = hl.read_matrix_table(get_ukb_pheno_mt_path('full'))
        meta_ht = hl.import_table(get_ukb_meta_pop_tsv_path(), impute=True, key='s')
        mt = mt.annotate_rows(**meta_ht[mt.row_key])

        ht = mt.group_rows_by('pop').aggregate(
            n_cases_by_pop=hl.cond(mt.data_type == 'continuous',
                                   hl.agg.count_where(hl.is_defined(mt.value)),
                                   hl.int64(hl.agg.sum(mt.value)))
        ).entries()
        ht.export(get_phenotype_summary_tsv_path('full'))

    if args.prepare_genotype_data:
        load_all_mfi_data().write(ukb_imputed_info_ht_path, args.overwrite)

    if args.genotype_summary:
        # variants = hl.read_table(ukb_imputed_info_ht_path)
        # print(variants.count())
        # variants = variants.filter(variants.info > 0.8)
        # print(variants.count())
        # meta_ht = hl.import_table(get_ukb_meta_pop_tsv_path(), impute=True, types={'s': hl.tstr}, key='s')
        # mt = get_ukb_imputed_data('all', variant_list=variants)
        # mt = mt.annotate_cols(**meta_ht[mt.col_key])
        # ht = mt.annotate_rows(af=hl.agg.group_by(mt.pop, hl.agg.mean(mt.dosage)),
        #                       an=hl.agg.group_by(mt.pop, hl.agg.count_where(hl.is_defined(mt.dosage)))).rows()
        ht = hl.utils.range_table(1)
        ht = ht.checkpoint(ukb_af_ht_path, args.overwrite, _read_if_exists=not args.overwrite)
        print(ht.aggregate(hl.struct(
            # hist=hl.agg.hist(hl.sum(ht.an.values()), 0, total_samples, 10),  # No missing data
            # fraction_missingness=hl.agg.fraction(hl.sum(ht.an.values()) < total_samples),
            number_sites_above_001=hl.agg.count_where(hl.any(lambda x: x / 2 > 0.001, ht.af.values())),
            number_sites_above_005=hl.agg.count_where(hl.any(lambda x: x / 2 > 0.005, ht.af.values())),
            number_sites_above_01=hl.agg.count_where(hl.any(lambda x: x / 2 > 0.01, ht.af.values())),
            number_sites_above_05=hl.agg.count_where(hl.any(lambda x: x / 2 > 0.05, ht.af.values())),
            number_sites_above_10=hl.agg.count_where(hl.any(lambda x: x / 2 > 0.1, ht.af.values()))
        )))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--pheno_summary', action='store_true')
    parser.add_argument('--prepare_genotype_data', action='store_true')
    parser.add_argument('--genotype_summary', action='store_true')
    parser.add_argument('--overwrite', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main)
    else:
        main()