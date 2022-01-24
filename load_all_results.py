#!/usr/bin/env python3

__author__ = 'konradk'

import sys
import argparse
from pprint import pprint
from datetime import date
from collections import Counter
from tqdm import tqdm
from gnomad.utils import slack
from ukbb_pan_ancestry import *


def get_all_valid_variant_results_ht_paths(pop):
    results_dir = f'{bucket}/results/result/{pop}'
    all_phenos_dir = hl.hadoop_ls(results_dir)
    all_variant_outputs = get_files_in_parent_directory(all_phenos_dir)
    return [x for x in all_variant_outputs if 'prescriptions--' not in x]


def patch_mt_keys(mt, patch_more: bool = False):
    mt = mt.key_cols_by(**{x: hl.case(missing_false=True)
                        .when((mt.phenocode == 'whr') & (x == 'modifier'), 'irnt')  # Patch for accidental coding of irnt-whr as "whr" for SAIGE run
                        .when((mt.trait_type == 'biomarkers') & (x == 'coding'), NULL_STR_KEY)  # Patch for accidental coding of biomarkers as their phenocode for SAIGE run
                        .when((mt.trait_type == 'biomarkers') & (x == 'modifier'), 'irnt')  # Patch for accidental coding of biomarkers as "" for SAIGE run
                        .when((mt.phenocode == '5097') & (x == 'modifier'), 'irnt')  # Patch for accidental coding of 5097 as "" for SAIGE run
                        .when(patch_more & (mt.trait_type == 'continuous') & (x == 'modifier') &
                              (mt.modifier == mt.phenocode), '')  # Patch for accidental coding of continuous as their phenocode for SAIGE run
                        .default(mt[x])
                           for x in PHENO_KEY_FIELDS})
    return mt


def apply_qc(mt, case_ac_threshold: int = 3, overall_mac_threshold: int = 20, min_case_count: int = 50):
    mt = mt.filter_cols(mt.n_cases >= min_case_count)
    ac_cases = mt.n_cases * 2 * mt['AF.Cases']
    an_controls = mt.n_controls * 2 * mt['AF.Controls']

    maf_total = 0.5 - hl.abs(0.5 - mt['AF_Allele2'])
    an_total = (mt.n_cases + hl.or_else(mt.n_controls, 0)) * 2
    mac_total = maf_total * an_total

    return mt.annotate_entries(
        low_confidence=hl.case(missing_false=True)
            .when(ac_cases <= case_ac_threshold, True)
            .when(an_controls <= case_ac_threshold, True)
            .when(mac_total <= overall_mac_threshold, True)
            .default(False)
    )


def check_and_annotate_with_dict(mt, input_dict, dict_key_from_mt, axis='cols'):
    direction = mt.col if axis == 'cols' else mt.row
    annotation = {}
    for key, value in list(input_dict.values()[0].items()):
        if key in list(direction):
            annotation[key] = hl.case().when(
                hl.is_defined(mt[key]), mt[key]).when(
                hl.is_defined(input_dict.get(dict_key_from_mt)), input_dict.get(dict_key_from_mt)[key]).default(
                hl.null(value.dtype)
            )
        else:
            annotation[key] = hl.if_else(hl.is_defined(input_dict.get(dict_key_from_mt)),
                                         input_dict.get(dict_key_from_mt)[key],
                                         hl.null(value.dtype))
    if axis == 'cols':
        mt = mt.annotate_cols(**annotation)
    else:
        mt = mt.annotate_rows(**annotation)
    return mt


def generate_sumstats_mt(all_variant_outputs, heritability_dict, pheno_dict, temp_dir, inner_mode = '_read_if_exists',
                         checkpoint: bool = False):
    row_keys = ['locus', 'alleles', 'gene', 'annotation']
    col_keys = PHENO_KEY_FIELDS

    all_hts = [unify_saige_ht_schema(hl.read_table(x), patch_case_control_count=x) for x in tqdm(all_variant_outputs)]
    mt = join_pheno_hts_to_mt(all_hts, row_keys, col_keys, temp_dir=temp_dir,
                              inner_mode=inner_mode, repartition_final=20000)

    if checkpoint:
        mt = mt.checkpoint(f'{temp_dir}/staging.mt', **{inner_mode: True})
    mt = patch_mt_keys(mt)
    key = mt.col_key.annotate(phenocode=format_pheno_dir(mt.phenocode))
    mt = check_and_annotate_with_dict(mt, pheno_dict, key)
    if mt.inv_normalized.dtype == hl.tstr:
        mt = mt.annotate_cols(inv_normalized=hl.bool(mt.inv_normalized))
    key = mt.col_key.annotate(phenocode=format_pheno_dir(mt.phenocode))
    mt = check_and_annotate_with_dict(mt, heritability_dict, key)
    mt = mt.filter_cols(mt.phenocode != "").drop('Is.SPA.converge', 'AC_Allele2')
    mt = mt.key_rows_by('locus', 'alleles')
    return mt


def write_full_mt(overwrite):
    mts = []
    for pop in POPS:
        mt = hl.read_matrix_table(get_variant_results_path(pop, 'mt')).annotate_cols(pop=pop)
        mt = mt.annotate_cols(_logged=hl.agg.any(mt.Pvalue < 0))
        mt = mt.annotate_entries(Pvalue=hl.if_else(mt._logged, mt.Pvalue, hl.log(mt.Pvalue))).drop('_logged')

        mt = mt.filter_cols((mt.coding != 'zekavat_20200409') & ~mt.phenocode.contains('covid'))
        mt = apply_qc(mt)
        mt = patch_mt_keys(mt)
        mt = reannotate_cols(mt, pop)
        mt = patch_mt_keys(mt, patch_more=True)
        mt = mt.filter_cols(hl.is_defined(mt.heritability) &
                            (hl.is_defined(mt.description) | (mt.trait_type == 'prescriptions')))
        mt = re_colkey_mt(mt)
        mt = mt.select_cols(pheno_data=mt.col_value)
        mt = mt.select_entries(summary_stats=mt.entry)
        # Fix for doubled phenos (new runs)
        mt = mt.collect_cols_by_key()
        def get_index(mt):
            return hl.sorted(hl.enumerate(mt.pheno_data), key=lambda x: x[1].saige_version.split('_')[1], reverse=True)[0][0]
        mt = mt.select_entries(summary_stats=mt.summary_stats[get_index(mt)])
        mt = mt.select_cols(pheno_data=mt.pheno_data[get_index(mt)])
        if pop == 'EUR':
            x = mt.filter_cols(mt.phenocode == '30610')
            print(x.aggregate_entries(hl.agg.counter(x.summary_stats.SE == 0)))
        mts.append(mt)

    full_mt = mts[0]
    for mt in mts[1:]:
        full_mt = full_mt.union_cols(mt, row_join_type='outer')
    full_mt = full_mt.collect_cols_by_key()
    full_mt = full_mt.annotate_cols(
        pheno_data=full_mt.pheno_data.map(lambda x: x.drop(*PHENO_COLUMN_FIELDS)),
        **{x: full_mt.pheno_data[x][0] for x in PHENO_DESCRIPTION_FIELDS},
        **{f'n_cases_full_cohort_{sex}': full_mt.pheno_data[f'n_cases_{sex}'][0]
           for sex in ('both_sexes', 'females', 'males')}
    )
    full_mt.write(get_variant_results_path('full'), overwrite)
    print('Pops per pheno:')
    pprint(dict(Counter(full_mt.aggregate_cols(hl.agg.counter(hl.len(full_mt.pheno_data))))))



def reannotate_cols(mt, pop):
    heritability_dict = get_heritability_dict(pop)
    pheno_dict = get_pheno_dict()
    key = get_modified_key(mt)
    mt = check_and_annotate_with_dict(mt, pheno_dict, key)
    key = get_modified_key(mt)
    return check_and_annotate_with_dict(mt, heritability_dict, key)


def recode_single_pheno_struct_to_legacy_path(pheno_struct):
    trait_type = pheno_struct.trait_type
    phenocode = pheno_struct.phenocode
    if trait_type == 'icd10':
        trait_type = 'icd_all'
        coding = 'icd10'
    elif trait_type == 'phecode':
        coding = pheno_struct.pheno_sex
    elif trait_type == 'biomarkers':
        coding = pheno_struct.phenocode
    else:
        if phenocode == 'whr':
            coding = 'whr'
        else:
            coding = pheno_struct.coding if pheno_struct.coding else pheno_struct.modifier
    return f'{trait_type}-{format_pheno_dir(phenocode)}-{coding}'


def main(args):
    hl.init(default_reference='GRCh37', log='/combine_results.log', branching_factor=8)

    inner_mode = 'overwrite' if args.overwrite else '_read_if_exists'
    pops = args.pops.split(',') if args.pops else POPS

    if args.run_basic_load:
        for pop in pops:
            all_variant_outputs = get_all_valid_variant_results_ht_paths(pop)
            pheno_dict = get_pheno_dict()
            heritability_dict = get_heritability_dict(pop)

            mt = generate_sumstats_mt(all_variant_outputs, heritability_dict, pheno_dict, f'{temp_bucket}/{pop}/variant', inner_mode)
            mt.write(get_variant_results_path(pop, 'mt'), overwrite=args.overwrite)

    if args.run_additional_load:
        today = date.today().strftime("%y%m%d")
        for pop in pops:
            all_variant_outputs = get_all_valid_variant_results_ht_paths(pop)
            pheno_dict = get_pheno_dict()
            heritability_dict = get_heritability_dict(pop)

            all_keys = hl.read_matrix_table(get_variant_results_path(pop, 'mt')).col_key.collect()
            loaded_phenos = {
                recode_single_pheno_struct_to_legacy_path(x) for x in all_keys
            }
            loaded_phenos |= {f'{x.trait_type}-{format_pheno_dir(x.phenocode)}-{x.pheno_sex}-{x.coding}-{x.modifier}' for x in all_keys}

            # Some phenotypes were manually modified later but had already made it in the release
            loaded_phenos |= {
                'categorical-41229-5N8', 'categorical-20004-1503', 'categorical-20004-1563',  # CSA
                'continuous-5097-irnt', 'continuous-104550-104550',  # EAS
            }

            def _matches_any_pheno(pheno_path, phenos_to_match):
                return any(x for x in phenos_to_match if f'/{x}/variant_results.ht' in pheno_path)

            if args.force_reload:
                pheno_matches = set(args.force_reload.split(','))
                all_variant_outputs = [x for x in all_variant_outputs if _matches_any_pheno(x, pheno_matches)]
            else:
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
                continue
            if not len(all_variant_outputs):
                continue

            mt = generate_sumstats_mt(all_variant_outputs, heritability_dict, pheno_dict,
                                      f'{temp_bucket}/{pop}/variant_{today}', inner_mode, checkpoint=True)

            original_mt = hl.read_matrix_table(get_variant_results_path(pop, 'mt'))
            original_mt = original_mt.checkpoint(f'{temp_bucket}/{pop}/variant_before_{today}.mt', overwrite=True)
            if args.force_reload:
                original = original_mt.count_cols()
                original_mt = original_mt.filter_cols(
                    hl.literal(pheno_matches).contains(hl.delimit(list(original_mt.col_key.values()), '-')), keep=False)
                print(f'\n\nGoing from {original} to {original_mt.count_cols()}...\n\n')
            mt = re_colkey_mt(mt)
            original_mt = re_colkey_mt(original_mt)
            mt = original_mt.union_cols(mt, row_join_type='outer')
            mt.write(get_variant_results_path(pop, 'mt'), overwrite=args.overwrite)

    if args.run_combine_load:
        write_full_mt(args.overwrite)


def re_colkey_mt(mt):
    mt = mt.key_cols_by().select_cols(*PHENO_KEY_FIELDS, *PHENO_COLUMN_FIELDS, *PHENO_GWAS_FIELDS, 'pop')
    return mt.key_cols_by(*PHENO_KEY_FIELDS)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--run_basic_load', help='Overwrite everything', action='store_true')
    parser.add_argument('--run_additional_load', help='Overwrite everything', action='store_true')
    parser.add_argument('--run_combine_load', help='Overwrite everything', action='store_true')
    parser.add_argument('--dry_run', help='Overwrite everything', action='store_true')
    parser.add_argument('--load_only', help='Comma-separated list of trait_type-pheno-coding to run'
                                            '(e.g. continuous-50-irnt,icd_all-E10-icd10 )')
    parser.add_argument('--force_reload', help='Comma-separated list of trait_type-pheno-coding to force reload'
                                            '(e.g. continuous-50-irnt,icd_all-E10-icd10 )')
    parser.add_argument('--find_errors', help='Overwrite everything', action='store_true')
    parser.add_argument('--pops', help='comma-separated list')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        from slack_token_pkg.slack_creds import slack_token
        with slack.slack_notifications(slack_token, args.slack_channel):
            main(args)
