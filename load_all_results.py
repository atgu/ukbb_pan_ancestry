#!/usr/bin/env python3

__author__ = 'konradk'

from ukb_common import *
from ukbb_pan_ancestry import *

bucket = 'gs://ukb-diverse-pops'
temp_bucket = 'gs://ukbb-diverse-temp-30day'


def main(args):
    hl.init(default_reference='GRCh38')

    POPS.remove('EUR')
    hts = []
    mts = []
    for pop in POPS:
        results_dir = f'{bucket}/results/result/{pop}'
        all_phenos = hl.hadoop_ls(results_dir)
        all_variant_outputs = []
        all_variant_mt_outputs = []
        for directory in all_phenos:
            variant_results_ht_path = f'{directory["path"]}/variant_results.ht'
            variant_results_mt_path = f'{directory["path"]}/variant_results.mt'
            if hl.hadoop_exists(f'{variant_results_mt_path}/_SUCCESS'):
                all_variant_outputs.append(variant_results_ht_path)
                all_variant_mt_outputs.append(variant_results_mt_path)
        if args.dry_run:
            print(f'Found {len(all_variant_mt_outputs)} inputs.')
            sys.exit(0)

        pheno_ht = hl.read_matrix_table(get_ukb_pheno_mt_path()).cols()
        pheno_dict = hl.dict(pheno_ht.aggregate(hl.agg.collect(
            (hl.struct(pheno=pheno_ht.pheno, coding=pheno_ht.coding, trait_type=pheno_ht.data_type), pheno_ht.row_value)), _localize=False))
        all_hts = list(map(lambda x: hl.read_table(x), all_variant_outputs))
        ht = all_hts[0].union(*all_hts[1:], unify=True)
        ht = ht.annotate(**pheno_dict[hl.struct(pheno=ht.pheno, coding=ht.coding, trait_type=ht.trait_type)])
        ht.write(get_variant_results_path(pop), overwrite=args.overwrite)

        raw_mts = list(map(lambda x: hl.read_matrix_table(x), all_variant_mt_outputs))
        all_mts = []
        for mt in raw_mts:
            if 'AF.Cases' not in list(mt.entry):
                mt = mt.select_entries('AC_Allele2', 'AF_Allele2', 'imputationInfo', 'N', 'BETA', 'SE', 'Tstat',
                                       **{'p.value.NA': hl.null(hl.tfloat64), 'Is.SPA.converge': hl.null(hl.tint32),
                                          'varT': mt.varT, 'varTstar': mt.varTstar, 'AF.Cases': hl.null(hl.tfloat64),
                                          'AF.Controls': hl.null(hl.tfloat64), 'Pvalue': mt.Pvalue})
            all_mts.append(mt)
        mt = union_mts_by_tree(all_mts, temp_bucket + '/variant')
        mt = mt.annotate_cols(**pheno_dict[mt.col_key])
        mt = mt.checkpoint(get_variant_results_path(pop, 'mt').replace(bucket, temp_bucket),
                           overwrite=args.overwrite, _read_if_exists=not args.overwrite)
        mt.repartition(5000, shuffle=False).write(get_variant_results_path(pop, 'mt'), overwrite=args.overwrite)

        mt = hl.read_matrix_table(get_variant_results_path(pop, 'mt'))
        mt = mt.annotate_cols(lambda_gc=hl.methods.statgen._lambda_gc_agg(mt.Pvalue))
        ht = mt.cols().annotate(pop=pop)
        hts.append(ht)
        ht.export(get_variant_results_path(pop, 'lambdas.txt.bgz'))

        mt = hl.read_matrix_table(get_variant_results_path(pop, 'mt')).annotate_cols(pop=pop)
        mt = mt.select_cols(pheno_data=mt.col_value)
        mt = mt.select_entries(summary_stats=mt.entry)
        mts.append(mt)

    full_ht = hts[0].union(*hts[1:])
    full_ht.write(f'{bucket}/combined_results/all_lambdas.ht')
    full_ht.group_by('pop').aggregate(lambda_gc=hl.agg.stats(full_ht.lambda_gc).drop('sum')).show(width=200)
    full_ht.group_by('pheno', 'meaning').aggregate(lambda_gc=hl.agg.stats(full_ht.lambda_gc).drop('sum')).show(20, width=200)

    full_mt = mts[0]
    for mt in mts[1:]:
        full_mt = full_mt.union_cols(mt, row_join_type='outer')
    full_mt = full_mt.collect_cols_by_key().annotate_globals(pops=POPS)
    full_mt.write(f'{bucket}/combined_results/all_sumstats.mt')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--dry_run', help='Overwrite everything', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)