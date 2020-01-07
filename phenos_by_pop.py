from gnomad_hail import *
from ukbb_pan_ancestry import *
from ukb_common import *

temp_bucket = 'gs://ukbb-diverse-temp-30day'


def main():
    data_types = ('categorical', 'continuous', 'icd')
    hl.init(default_reference='GRCh37', log='/pheno_agg.log')

    for data_type in data_types:
        mt = hl.read_matrix_table(get_ukb_pheno_mt_path(data_type))
        meta_ht = hl.import_table(get_ukb_meta_pop_tsv_path(), impute=True, key='s')
        mt = mt.annotate_rows(**meta_ht[mt.row_key])

        if data_type == 'continuous':
            mt = mt.filter_cols(mt.coding != 'raw')

        extra_fields = compute_n_cases(mt, data_type)
        ht = mt.group_rows_by('pop').aggregate(**extra_fields).entries()
        ht.export(get_phenotype_summary_tsv_path(data_type))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main)
    else:
        main()