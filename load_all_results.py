#!/usr/bin/env python3

__author__ = 'konradk'

import sys
import argparse
from pprint import pprint
from datetime import date
from collections import defaultdict, Counter
from ukbb_pan_ancestry import *
from tqdm import tqdm


def get_all_valid_variant_results_ht_paths(pop):
    results_dir = f'{bucket}/results/result/{pop}'
    all_phenos_dir = hl.hadoop_ls(results_dir)
    all_variant_outputs = get_files_in_parent_directory(all_phenos_dir)
    return [x for x in all_variant_outputs if 'prescriptions--' not in x]


def get_heritability_dict(pop):
    heritability_ht = hl.import_table(get_heritability_txt_path(), impute=True, key=PHENO_KEY_FIELDS)
    heritability_ht = heritability_ht.filter(heritability_ht.pop == pop)
    heritability_dict = create_broadcast_dict(heritability_ht.key)
    return heritability_dict


def generate_lambda_ht(pop, k: int = 75):
    _intervals = get_n_even_intervals(5000)
    mt = hl.read_matrix_table(get_variant_results_path(pop, 'mt'))#, _intervals=_intervals)
    ac_cutoffs = list(range(0, 6)) + [10, 20, 50, 100]
    af_cutoffs = sorted([0] + [y * 10 ** x for y in (1, 2, 5) for x in range(-4, 0)] + [0.99])

    af_cases = mt['AF.Cases']
    ac_cases = af_cases * mt.n_cases * 2

    af_total = mt['AF_Allele2']
    ac_total = af_total * 2 * (mt.n_cases + mt.n_controls)
    p_value_field = mt.Pvalue

    sig_threshold = 5e-8

    mt = mt.annotate_cols(sumstats_qc=hl.struct(**{
        f'{metric}_{breakdown}_{flavor}': [hl.agg.filter(breakdown_dict[flavor] >= cutoff, agg) for cutoff in cutoffs]
        for flavor, cutoffs in (('ac', ac_cutoffs), ('af', af_cutoffs))
        for breakdown, breakdown_dict in (('by_case', {'ac': ac_cases}),  #, 'af': af_cases}),
                                          ('by', {# 'ac': ac_total,
                                                  'af': af_total})) if flavor in breakdown_dict
        for metric, agg in (
            ('lambda_gc', hl.methods.statgen._lambda_gc_agg(p_value_field, k=k)),
            ('n_variants', hl.agg.count()),
            ('n_sig', hl.agg.count_where(p_value_field < sig_threshold))
        )
    })).annotate_globals(ac_cutoffs=ac_cutoffs, af_cutoffs=af_cutoffs)

    ht = mt.cols()
    ht = ht.key_by(pop=pop, *PHENO_KEY_FIELDS)
    return ht


def generate_sumstats_mt(all_variant_outputs, heritability_dict, pheno_dict, temp_dir, inner_mode = '_read_if_exists'):
    row_keys = ['locus', 'alleles', 'gene', 'annotation']
    col_keys = PHENO_KEY_FIELDS

    all_hts = [unify_saige_ht_schema(hl.read_table(x), patch_case_control_count=x) for x in tqdm(all_variant_outputs)]
    mt = join_pheno_hts_to_mt(all_hts, row_keys, col_keys, temp_dir=temp_dir,
                              inner_mode=inner_mode, repartition_final=20000)
    # Patch for accidental coding of irnt-whr as "whr" for SAIGE run
    mt = mt.key_cols_by(**{x: hl.if_else((mt.phenocode == 'whr') & (x == 'modifier'), 'irnt', mt[x])
                           for x in PHENO_KEY_FIELDS})
    key = mt.col_key.annotate(phenocode=format_pheno_dir(mt.phenocode))
    mt = mt.annotate_cols(**pheno_dict.get(key), **heritability_dict.get(key))
    mt = mt.filter_cols(mt.phenocode != "").drop('varT', 'varTstar', 'Is.SPA.converge',
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
    description_fields = ('description', 'description_more', 'coding_description', 'category')
    full_mt = full_mt.annotate_cols(
        pheno_data=full_mt.pheno_data.map(lambda x: x.drop(*description_fields,
                                                           'n_cases_both_sexes', 'n_cases_males', 'n_cases_females')),
        **{x: full_mt.pheno_data[x][0] for x in description_fields},
        **{f'n_cases_full_cohort_{sex}': full_mt.pheno_data[f'n_cases_{sex}'][0]
           for sex in ('both_sexes', 'females', 'males')}
    )
    full_mt.write(get_variant_results_path('full'), overwrite)


def generate_and_write_lambdas(overwrite):
    lambda_hts = []
    for pop in POPS:
        ht = generate_lambda_ht(pop)
        ht = ht.checkpoint(get_variant_results_path(pop, 'lambdas.ht'), overwrite=overwrite, _read_if_exists=not overwrite)
        explode_lambda_ht(ht).export(get_variant_results_path(pop, 'lambdas_by_ac.txt.bgz'))
        explode_lambda_ht(ht, 'af').export(get_variant_results_path(pop, 'lambdas_by_af.txt.bgz'))
        lambda_hts.append(ht)

    full_ht = lambda_hts[0].union(*lambda_hts[1:])
    pheno_ht = hl.read_table(get_phenotype_summary_path('full')).select('n_cases_by_pop', 'stats')
    full_ht = full_ht.annotate(**pheno_ht[full_ht.key])
    full_ht = full_ht.checkpoint(f'{bucket}/combined_results/all_lambdas.ht', overwrite=overwrite, _read_if_exists=not overwrite)
    # full_ht.group_by('pop').aggregate(lambda_gc=hl.agg.stats(full_ht.lambda_gc).drop('sum')).show(width=200)
    # full_ht.group_by('pheno', 'meaning').aggregate(lambda_gc=hl.agg.stats(full_ht.lambda_gc).drop('sum')).show(20, width=200)
    # full_ht.group_by('pop').aggregate(n_samples=hl.agg.max(full_ht.n_cases_by_pop), n_phenos=hl.agg.count()).show()
    print('Pops per pheno:')
    pprint(dict(Counter(full_ht.aggregate(hl.agg.group_by(full_ht.key.drop('pop'), hl.agg.count()).values()))))
    full_ht.flatten().export(f'{bucket}/combined_results/all_lambdas.txt.bgz')


def explode_lambda_ht(ht, by='ac'):
    ac_ht = ht.annotate(sumstats_qc=ht.sumstats_qc.select(*[x for x in ht.sumstats_qc.keys() if f'_{by}' in x]))
    ac_ht = ac_ht.annotate(index_ac=hl.zip_with_index(ac_ht[f'{by}_cutoffs'])).explode('index_ac')
    ac_ht = ac_ht.transmute(**{by: ac_ht.index_ac[1]},
                            **{x: ac_ht.sumstats_qc[x][ac_ht.index_ac[0]] for x in ac_ht.sumstats_qc})
    return ac_ht


def main(args):
    hl.init(default_reference='GRCh37', log='/combine_results.log')

    inner_mode = 'overwrite' if args.overwrite else '_read_if_exists'
    pops = args.pops.split(',') if args.pops else POPS

    if args.run_basic_load:
        for pop in pops:
            all_variant_outputs = get_all_valid_variant_results_ht_paths(pop)
            pheno_dict = get_pheno_dict()
            heritability_dict = get_heritability_dict(pop)

            mt = generate_sumstats_mt(all_variant_outputs, heritability_dict, pheno_dict, f'{temp_bucket}/{pop}/variant', inner_mode)
            mt.write(get_variant_results_path(pop, 'mt'), overwrite=args.overwrite)
        if pops != POPS: sys.exit(0)

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
                    print('WARNING: Empty string in pheno_matches. Might reload more than expected')
                all_variant_outputs = [x for x in all_variant_outputs if _matches_any_pheno(x, pheno_matches)]

            print(f'Loading {len(all_variant_outputs)} additional HTs...')
            if len(all_variant_outputs) < 20:
                print(all_variant_outputs)
            if args.dry_run:
                sys.exit(0)

            mt = generate_sumstats_mt(all_variant_outputs, heritability_dict, pheno_dict,
                                      f'{temp_bucket}/{pop}/variant_{today}', inner_mode)

            original_mt = hl.read_matrix_table(get_variant_results_path(pop, 'mt'))
            original_mt = original_mt.checkpoint(f'{temp_bucket}/{pop}/variant_before_{today}.mt', _read_if_exists=True)
            mt = original_mt.union_cols(mt, row_join_type='outer')
            mt.write(get_variant_results_path(pop, 'mt'), overwrite=args.overwrite)

    # _intervals = get_n_even_intervals(5000)
    # print(_intervals[0])

    if args.test_approx_median_k:
        # test_approx_median_k()
        # ht = hl.read_matrix_table(get_variant_results_path('AFR', 'mt')).rows()
        # ht = ht.collect_by_key()
        # ht = ht.filter((hl.len(ht.values) > 1))
        # pprint(ht.aggregate(hl.agg.counter(hl.delimit(hl.sorted(ht.values.annotation)))))
        ht = hl.read_matrix_table(get_variant_results_path('AFR', 'mt')).rows()
        print(ht._force_count())
        ht = hl.read_matrix_table(get_variant_results_path('AFR', 'mt'), _intervals=_intervals).rows()
        print(ht._force_count())

    if args.run_basic_load or args.run_additional_load or args.run_combine_load:
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


def test_approx_median_k():
    ht = hl.read_table('gs://ukb-diverse-pops/results/result/AFR/icd_all-I48-icd10/variant_results.ht')
    ks = [x * 10 for x in range(1, 11)]
    ac_cutoffs = list(range(0, 6)) + [10, 20, 50, 100]
    af_cutoffs = sorted([0] + [y * 10 ** x for y in (1, 2, 5) for x in range(-4, 0)] + [0.99])
    af_cases = ht['AF.Cases']
    ac_cases = af_cases * ht.n_cases * 2
    af_total = ht['AF_Allele2']
    ac_total = af_total * 2 * (ht.n_cases + ht.n_controls)
    p_value_field = ht.Pvalue
    sig_threshold = 5e-8
    test = ht.group_by(_x=True).aggregate(
        sumstats_qc=hl.struct(**{
            f'{metric}_{breakdown}_{flavor}_{k}': [hl.agg.filter(breakdown_dict[flavor] >= cutoff, agg) for cutoff in
                                                   cutoffs]
            for flavor, cutoffs in (('ac', ac_cutoffs), ('af', af_cutoffs))
            for breakdown, breakdown_dict in (('by_case', {'ac': ac_cases}),  # , 'af': af_cases}),
                                              ('by', {  # 'ac': ac_total,
                                                  'af': af_total})) if flavor in breakdown_dict
            for k in ks
            for metric, agg in (
                ('lambda_gc', hl.methods.statgen._lambda_gc_agg(p_value_field, k=k)),
                ('n_variants', hl.agg.count()),
                ('n_sig', hl.agg.count_where(p_value_field < sig_threshold))
            )
        })
    ).annotate_globals(ac_cutoffs=ac_cutoffs, af_cutoffs=af_cutoffs)
    test = test.checkpoint(f'{temp_bucket}/test.ht', overwrite=True)
    explode_lambda_ht(test).export(f'{temp_bucket}/test_ac.txt.bgz')
    explode_lambda_ht(test, 'af').export(f'{temp_bucket}/test_af.txt.bgz')


def get_pheno_dict():
    pheno_ht = hl.read_matrix_table(get_ukb_pheno_mt_path()).cols()
    pheno_dict = create_broadcast_dict(pheno_ht.key)
    return pheno_dict


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--run_basic_load', help='Overwrite everything', action='store_true')
    parser.add_argument('--run_additional_load', help='Overwrite everything', action='store_true')
    parser.add_argument('--run_combine_load', help='Overwrite everything', action='store_true')
    parser.add_argument('--test_approx_median_k', help='Overwrite everything', action='store_true')
    parser.add_argument('--dry_run', help='Overwrite everything', action='store_true')
    parser.add_argument('--load_only', help='Comma-separated list of trait_type-pheno-coding to run'
                                            '(e.g. continuous-50-irnt,icd_all-E10-icd10 )')
    parser.add_argument('--find_errors', help='Overwrite everything', action='store_true')
    parser.add_argument('--pops', help='comma-separated list')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        from slack_token_pkg.slack_creds import slack_token
        with slack.slack_notifications(slack_token, args.slack_channel):
            main(args)
