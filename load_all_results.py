#!/usr/bin/env python3

__author__ = 'konradk'

import sys
import argparse
from pprint import pprint
from datetime import date
from ukb_common import *
from ukbb_pan_ancestry import *


def get_all_valid_variant_results_ht_paths(pop):
    results_dir = f'{bucket}/results/result/{pop}'
    all_phenos_dir = hl.hadoop_ls(results_dir)
    all_variant_outputs = get_files_in_parent_directory(all_phenos_dir)
    return [x for x in all_variant_outputs if 'prescriptions--' not in x]


def get_heritability_dict(pop):
    heritability_ht = hl.import_table(get_heritability_txt_path(), impute=True)
    heritability_ht = heritability_ht.filter(heritability_ht.pop == pop)
    heritability_dict = create_broadcast_dict(
        hl.struct(pheno=heritability_ht.pheno, coding=heritability_ht.coding, trait_type=heritability_ht.trait),
        hl.struct(heritability=heritability_ht.heritability, inv_normalized=heritability_ht.inv_normalized))
    return heritability_dict


def generate_lambda_ht(pop):
    mt = hl.read_matrix_table(get_variant_results_path(pop, 'mt'))
    mt = mt.annotate_cols(lambda_gc=hl.methods.statgen._lambda_gc_agg(mt.Pvalue),
                          n_vars=hl.agg.count_where(hl.is_defined(mt.Pvalue)))
    ht = mt.cols()
    ht = ht.key_by(pop=pop, pheno=ht.pheno, coding=ht.coding, trait_type=ht.trait_type)
    return ht


def generate_sumstats_mt(all_variant_outputs, heritability_dict, pheno_dict, temp_dir, inner_mode = '_read_if_exists'):
    row_keys = ['locus', 'alleles', 'gene', 'annotation']
    col_keys = ['pheno', 'coding', 'trait_type']

    all_hts = list(map(lambda x: unify_saige_ht_schema(hl.read_table(x), patch_case_control_count=x), all_variant_outputs))
    mt = join_pheno_hts_to_mt(all_hts, row_keys, col_keys, pheno_dict, temp_dir,
                              inner_mode=inner_mode, repartition_final=20000)
    entry = heritability_dict.get(mt.col_key)
    mt = mt.annotate_cols(saige_heritability=entry.heritability)  #, saige_inv_normalized=entry.inv_normalized)
    mt = mt.filter_cols(mt.pheno != "").drop('varT', 'varTstar', 'Is.SPA.converge',
                                             'AC_Allele2', 'N', 'Tstat')
    mt = mt.key_rows_by('locus', 'alleles')
    return mt


def write_full_mt(overwrite):
    mts = []
    for pop in POPS:
        mt = hl.read_matrix_table(get_variant_results_path(pop, 'mt')).annotate_cols(pop=pop)
        mt = mt.select_cols(pheno_data=mt.col_value)
        mt = mt.select_entries(summary_stats=mt.entry)
        mts.append(mt)

    full_mt = mts[0]
    for mt in mts[1:]:
        full_mt = full_mt.union_cols(mt, row_join_type='outer')
    full_mt = full_mt.collect_cols_by_key()
    full_mt = full_mt.annotate_cols(
        pheno_data=full_mt.pheno_data.map(lambda x: x.drop('data_type', 'meaning', 'path',
                                                           'n_cases_both_sexes', 'n_cases_males', 'n_cases_females')),
        meaning=full_mt.pheno_data.meaning[0],
        description=full_mt.pheno_data.path[0],
        **{f'n_cases_full_cohort_{sex}': full_mt.pheno_data[f'n_cases_{sex}'][0]
           for sex in ('both_sexes', 'females', 'males')}
    )
    full_mt.write(get_variant_results_path('full'), overwrite)


def generate_and_write_lambdas(overwrite):
    lambda_hts = []
    for pop in POPS:
        ht = generate_lambda_ht(pop)
        ht = ht.checkpoint(get_variant_results_path(pop, 'lambdas.ht'), overwrite=overwrite, _read_if_exists=not overwrite)
        ht.export(get_variant_results_path(pop, 'lambdas.txt.bgz'))
        lambda_hts.append(ht)

    full_ht = lambda_hts[0].union(*lambda_hts[1:])
    pheno_ht = hl.read_table(get_phenotype_summary_path('full')).select('n_cases_by_pop', 'stats')
    full_ht = full_ht.annotate(**pheno_ht[full_ht.key])
    full_ht = full_ht.checkpoint(f'{bucket}/combined_results/all_lambdas.ht', overwrite=overwrite, _read_if_exists=not overwrite)
    full_ht.group_by('pop').aggregate(lambda_gc=hl.agg.stats(full_ht.lambda_gc).drop('sum')).show(width=200)
    full_ht.group_by('pheno', 'meaning').aggregate(lambda_gc=hl.agg.stats(full_ht.lambda_gc).drop('sum')).show(20, width=200)
    full_ht.group_by('pop').aggregate(n_samples=hl.agg.max(full_ht.n_cases_by_pop), n_phenos=hl.agg.count()).show()
    print('Pops per pheno:')
    pprint(dict(Counter(full_ht.aggregate(hl.agg.group_by(full_ht.key.drop('pop'), hl.agg.count()).values()))))
    full_ht.flatten().export(f'{bucket}/combined_results/all_lambdas.txt.bgz')


def main(args):
    hl.init(default_reference='GRCh37', log='/combine_results.log')

    inner_mode = 'overwrite' if args.overwrite else '_read_if_exists'

    if args.run_basic_load:
        for pop in POPS:
            all_variant_outputs = get_all_valid_variant_results_ht_paths(pop)
            pheno_dict = get_pheno_dict()
            heritability_dict = get_heritability_dict(pop)

            mt = generate_sumstats_mt(all_variant_outputs, heritability_dict, pheno_dict, f'{temp_bucket}/{pop}/variant', inner_mode)
            mt.write(get_variant_results_path(pop, 'mt'), overwrite=args.overwrite)

    if args.run_additional_load:
        today = date.today().strftime("%y%m%d")
        for pop in POPS:
            all_variant_outputs = get_all_valid_variant_results_ht_paths(pop)
            pheno_dict = get_pheno_dict()
            heritability_dict = get_heritability_dict(pop)

            loaded_phenos = set(
                [f'{x.trait_type}-{format_pheno_dir(x.pheno)}-{x.coding}' for x in hl.read_matrix_table(get_variant_results_path(pop, 'mt')).col_key.collect()]
            )

            def _matches_any_pheno(pheno_path, phenos_to_match):
                return any(x for x in phenos_to_match if f'/{x}/variant_results.ht' in pheno_path)

            all_variant_outputs = [x for x in all_variant_outputs if not _matches_any_pheno(x, loaded_phenos)]

            if args.load_only:
                pheno_matches = set(args.load_only.split(','))
                if '' in pheno_matches:
                    logger.warning('Empty string in pheno_matches. Might reload more than expected')
                all_variant_outputs = [x for x in all_variant_outputs if _matches_any_pheno(x, pheno_matches)]

            print(f'Loading {len(all_variant_outputs)} additional HTs...')
            if len(all_variant_outputs) < 20:
                print(all_variant_outputs)
            if args.dry_run:
                sys.exit(0)

            mt = generate_sumstats_mt(all_variant_outputs, heritability_dict, pheno_dict,
                                      f'{temp_bucket}/{pop}/variant_{today}', inner_mode)

            original_mt = hl.read_matrix_table(get_variant_results_path(pop, 'mt'))
            original_mt = original_mt.checkpoint(f'{temp_bucket}/{pop}/variant_before_{today}.mt', overwrite=args.overwrite)
            mt = original_mt.union_cols(mt)
            mt.write(get_variant_results_path(pop, 'mt'), overwrite=args.overwrite)


    if args.run_basic_load or args.run_additional_load:
        generate_and_write_lambdas(args.overwrite)
        write_full_mt(args.overwrite)


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


def get_pheno_dict():
    pheno_ht = hl.read_matrix_table(get_ukb_pheno_mt_path()).cols()
    pheno_dict = create_broadcast_dict(
        hl.struct(pheno=pheno_ht.pheno, coding=pheno_ht.coding, trait_type=pheno_ht.data_type),
        pheno_ht.row_value)
    return pheno_dict


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--run_basic_load', help='Overwrite everything', action='store_true')
    parser.add_argument('--run_additional_load', help='Overwrite everything', action='store_true')
    parser.add_argument('--dry_run', help='Overwrite everything', action='store_true')
    parser.add_argument('--load_only', help='Comma-separated list of trait_type-pheno-coding to run'
                                            '(e.g. continuous-50-irnt,icd_all-E10-icd10 )')
    parser.add_argument('--find_errors', help='Overwrite everything', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)