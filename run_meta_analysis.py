#!/usr/bin/env python
# coding: utf-8
import argparse
import hail as hl
from ukbb_pan_ancestry.resources.results import *
from ukbb_pan_ancestry.utils.results import *


def all_and_leave_one_out(x, pop_array, all_f=hl.sum, loo_f=lambda i, x: hl.sum(x) - hl.or_else(x[i], 0)):
    """
    Applies a function to an input array for all populations, and for each of leave-one-out populations.

    :param x: Input array
    :param pop_array: Population array
    :param all_f: Function for all populations. It takes the input array and returns a new value
    :param loo_f: Function for each of leave-one-out populations. It takes an index of leave-one-out
                  population and the input array, and returns an array of new values.
    ...
    :return: Array of new values for all populations and for each of leave-one-out populations.
    :rtype: ArrayExpression
    """
    arr = hl.array([all_f(x)])
    arr = arr.extend(hl.map(lambda i: loo_f(i, x), hl.range(hl.len(pop_array))))
    return hl.or_missing(hl.any(hl.is_defined, x), arr)


def main(args):
    hl.init()

    # Read in all sumstats
    mt = load_final_sumstats_mt(filter_phenos=True,
                                filter_variants=False,
                                filter_sumstats=True,
                                separate_columns_by_pop=False,
                                annotate_with_nearest_gene=False)

    # Annotate per-entry sample size
    def get_n(pheno_data, i):
        return pheno_data[i].n_cases + hl.or_else(pheno_data[i].n_controls, 0)

    mt = mt.annotate_entries(
        summary_stats=hl.map(lambda x: x[1].annotate(N=hl.or_missing(hl.is_defined(x[1]), get_n(mt.pheno_data, x[0]))),
                             hl.zip_with_index(mt.summary_stats)))

    # Exclude entries with low confidence flag.
    if not args.keep_low_confidence_variants:
        mt = mt.annotate_entries(summary_stats=hl.map(lambda x: hl.or_missing(~x.low_confidence, x), mt.summary_stats))

    # Run fixed-effect meta-analysis (all + leave-one-out)
    mt = mt.annotate_entries(unnorm_beta=mt.summary_stats.BETA / (mt.summary_stats.SE**2),
                             inv_se2=1 / (mt.summary_stats.SE**2))
    mt = mt.annotate_entries(sum_unnorm_beta=all_and_leave_one_out(mt.unnorm_beta, mt.pheno_data.pop),
                             sum_inv_se2=all_and_leave_one_out(mt.inv_se2, mt.pheno_data.pop))
    mt = mt.transmute_entries(META_BETA=mt.sum_unnorm_beta / mt.sum_inv_se2,
                              META_SE=hl.map(lambda x: hl.sqrt(1 / x), mt.sum_inv_se2))
    mt = mt.annotate_entries(META_Pvalue=hl.map(lambda x: 2 * hl.pnorm(x), -hl.abs(mt.META_BETA / mt.META_SE)))

    # Run heterogeneity test (Cochran's Q)
    mt = mt.annotate_entries(META_Q=hl.map(lambda x: hl.sum((mt.summary_stats.BETA - x)**2 * mt.inv_se2), mt.META_BETA),
                             variant_exists=hl.map(lambda x: ~hl.is_missing(x), mt.summary_stats.BETA))
    mt = mt.annotate_entries(META_N_pops=all_and_leave_one_out(mt.variant_exists, mt.pheno_data.pop))
    mt = mt.annotate_entries(META_Pvalue_het=hl.map(lambda i: hl.pchisqtail(mt.META_Q[i], mt.META_N_pops[i] - 1),
                                                    hl.range(hl.len(mt.META_Q))))

    # Add other annotations
    mt = mt.annotate_entries(ac_cases=hl.map(lambda x: x["AF.Cases"] * x.N, mt.summary_stats),
                             ac_controls=hl.map(lambda x: x["AF.Controls"] * x.N, mt.summary_stats),
                             META_AC_Allele2=all_and_leave_one_out(mt.summary_stats.AF_Allele2 * mt.summary_stats.N,
                                                                   mt.pheno_data.pop),
                             META_N=all_and_leave_one_out(mt.summary_stats.N, mt.pheno_data.pop))
    mt = mt.annotate_entries(META_AF_Allele2=mt.META_AC_Allele2 / mt.META_N,
                             META_AF_Cases=all_and_leave_one_out(mt.ac_cases, mt.pheno_data.pop) / mt.META_N,
                             META_AF_Controls=all_and_leave_one_out(mt.ac_controls, mt.pheno_data.pop) / mt.META_N)
    mt = mt.drop('unnorm_beta', 'inv_se2', 'variant_exists', 'ac_cases', 'ac_controls', 'summary_stats',
                 'META_AC_Allele2')

    # Format everything into array<struct>
    def is_finite_or_missing(x):
        return (hl.or_missing(hl.is_finite(x), x))

    meta_fields = ['BETA', 'SE', 'Pvalue', 'Q', 'Pvalue_het', 'N', 'N_pops', 'AF_Allele2', 'AF_Cases', 'AF_Controls']
    mt = mt.transmute_entries(meta_analysis=hl.map(
        lambda i: hl.struct(**{field: is_finite_or_missing(mt[f'META_{field}'][i]) for field in meta_fields}),
        hl.range(hl.len(mt.META_BETA))))

    col_fields = ['n_cases', 'n_controls']
    mt = mt.annotate_cols(
        **{field: all_and_leave_one_out(mt.pheno_data[field], mt.pheno_data.pop) for field in col_fields})
    col_fields += ['pop']
    mt = mt.annotate_cols(pop=all_and_leave_one_out(
        mt.pheno_data.pop,
        mt.pheno_data.pop,
        all_f=lambda x: x,
        loo_f=lambda i, x: hl.filter(lambda y: y != x[i], x),
    ))
    mt = mt.transmute_cols(meta_analysis_data=hl.map(
        lambda i: hl.struct(**{field: mt[field][i] for field in col_fields}), hl.range(hl.len(mt.pop))))

    mt.describe()
    mt.write(get_meta_analysis_results_path(), overwrite=args.overwrite)

    hl.copy_log('gs://ukb-diverse-pops/combined_results/meta_analysis.log')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--keep-low-confidence-variants',
                        help='Keep variants with low confidence flag for meta-analysis',
                        action='store_true')
    parser.add_argument('--overwrite', help='Overwrite data', action='store_true')
    args = parser.parse_args()

    main(args)
