#!/usr/bin/env python3

__author__ = 'konradk'

import argparse
from ukb_common import *
from ukbb_pan_ancestry import *

bucket = 'gs://ukb-diverse-pops'
temp_bucket = 'gs://ukbb-diverse-temp-30day'


def main(args):
    hl.init(default_reference='GRCh37')

    lambda_hts = []
    mts = []
    row_keys = ['locus', 'alleles', 'gene', 'annotation']
    col_keys = ['pheno', 'coding', 'trait_type']
    for pop in POPS:
        if pop == 'EUR': continue
        results_dir = f'{bucket}/results/result/{pop}'
        all_phenos_dir = hl.hadoop_ls(results_dir)
        all_variant_outputs = get_files_in_parent_directory(all_phenos_dir)

        pheno_ht = hl.read_matrix_table(get_ukb_pheno_mt_path()).cols()
        pheno_dict = hl.dict(pheno_ht.aggregate(hl.agg.collect(
            (hl.struct(pheno=pheno_ht.pheno, coding=pheno_ht.coding, trait_type=pheno_ht.data_type), pheno_ht.row_value)), _localize=False))

        all_hts = list(map(lambda x: unify_saige_ht_schema(hl.read_table(x)), all_variant_outputs))
        mt = join_pheno_hts_to_mt(all_hts, row_keys, col_keys, pheno_dict, f'{temp_bucket}/{pop}/variant',
                                  inner_mode='overwrite', repartition_final=20000)
        mt = pull_out_fields_from_entries(mt, ['AC_Allele2', 'AF_Allele2', 'imputationInfo', 'N']
                                          ).drop('varT', 'varTstar', 'Is.SPA.converge')
        mt = mt.filter_cols(mt.pheno != "")
        mt.key_rows_by('locus', 'alleles').write(get_variant_results_path(pop, 'mt'), overwrite=args.overwrite)

        mt = hl.read_matrix_table(get_variant_results_path(pop, 'mt'))
        mt = mt.annotate_cols(lambda_gc=hl.methods.statgen._lambda_gc_agg(mt.Pvalue),
                              n_vars=hl.agg.count_where(hl.is_defined(mt.Pvalue)))
        ht = mt.cols()
        ht = ht.key_by(pop=pop, pheno=ht.pheno, coding=ht.coding, trait_type=ht.trait_type)
        ht = ht.checkpoint(get_variant_results_path(pop, 'lambdas.ht'), overwrite=args.overwrite,
                           _read_if_exists=not args.overwrite)
        lambda_hts.append(ht)
        ht.export(get_variant_results_path(pop, 'lambdas.txt.bgz'))

        mt = hl.read_matrix_table(get_variant_results_path(pop, 'mt')).annotate_cols(pop=pop)
        mt = mt.select_cols(pheno_data=mt.col_value)
        mt = mt.select_entries(summary_stats=mt.entry)
        mts.append(mt)
    return

    full_ht = lambda_hts[0].union(*lambda_hts[1:])
    pheno_ht = hl.read_table(get_phenotype_summary_path('full')).select('n_cases_by_pop', 'stats')
    full_ht = full_ht.annotate(**pheno_ht[full_ht.key])
    full_ht = full_ht.checkpoint(f'{bucket}/combined_results/all_lambdas.ht', overwrite=args.overwrite, _read_if_exists=not args.overwrite)
    full_ht.group_by('pop').aggregate(lambda_gc=hl.agg.stats(full_ht.lambda_gc).drop('sum')).show(width=200)
    full_ht.group_by('pheno', 'meaning').aggregate(lambda_gc=hl.agg.stats(full_ht.lambda_gc).drop('sum')).show(20, width=200)

    full_mt = mts[0]
    for mt in mts[1:]:
        full_mt = full_mt.union_cols(mt, row_join_type='outer')
    full_mt = full_mt.collect_cols_by_key().annotate_globals(pops=POPS)
    # TODO: get proper case control counts in here
    full_mt.write(get_variant_results_path('full', 'mt'), args.overwrite)

    if args.find_errors:
        for pop in POPS:
            results_dir = f'{bucket}/results/result/{pop}'
            all_phenos = hl.hadoop_ls(results_dir)
            pheno_dirs = [x for x in all_phenos if x['is_dir']]
            all_errors = defaultdict(dict)
            for directory in pheno_dirs[1:]:
                _, pheno = directory['path'].rsplit('/', 1)
                all_files = hl.hadoop_ls(directory['path'])
                log_files = [x['path'] for x in all_files if x['path'].endswith('.log')]
                errors = hl.grep('[Ee]rror', log_files, show=False)
                for fname, errstrings in errors.items():
                    for errstring in errstrings:
                        if pheno not in all_errors[errstring]:
                            all_errors[errstring][pheno] = []
                        all_errors[errstring][pheno].append(fname)
            pprint(dict(all_errors).keys())


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--find_errors', help='Overwrite everything', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)