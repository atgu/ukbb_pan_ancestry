#!/usr/bin/env python
# coding: utf-8
import argparse
import hail as hl
from ukbb_pan_ancestry import *
from ukb_common import *


def all_and_leave_one_out(x, pop, all_f=hl.sum, loo_f=lambda i, x: hl.sum(x) - hl.or_else(x[i], 0)):
    arr = hl.array([all_f(x)])
    arr = arr.extend(hl.map(lambda i: loo_f(i, x), hl.range(hl.len(pop))))
    return arr


def main(args):
    hl.init(log='/run_meta_analysis.log')

    # Read in all sumstats
    mt = hl.read_matrix_table('gs://ukb-diverse-pops/combined_results/all_sumstats.mt')

    # Run meta-analysis (all + leave-one-out)
    mt = mt.annotate_entries(unnorm_beta=mt.summary_stats.BETA / (mt.summary_stats.SE**2),
                             inv_se2=1 / (mt.summary_stats.SE**2))

    mt = mt.transmute_entries(sum_unnorm_beta=all_and_leave_one_out(mt.unnorm_beta, mt.pheno_data.pop),
                              sum_inv_se2=all_and_leave_one_out(mt.inv_se2, mt.pheno_data.pop))

    mt = mt.transmute_entries(META_BETA=mt.sum_unnorm_beta / mt.sum_inv_se2,
                              META_SE=hl.map(lambda x: hl.sqrt(1 / x), mt.sum_inv_se2))

    mt = mt.annotate_entries(META_BETA=hl.map(lambda x: hl.or_missing(hl.is_finite(x), x), mt.META_BETA),
                             META_SE=hl.map(lambda x: hl.or_missing(hl.is_finite(x), x), mt.META_SE))

    mt = mt.annotate_entries(META_Pvalue=hl.map(lambda x: 2 * hl.pnorm(x), -hl.abs(mt.META_BETA / mt.META_SE)))

    # Add other annotations
    mt = mt.annotate_entries(variant_exists=hl.map(lambda x: ~hl.is_missing(x), mt.summary_stats.BETA),
                             af_cases=hl.map(lambda x: x["AF.Cases"] * x.N, mt.summary_stats),
                             af_controls=hl.map(lambda x: x["AF.Controls"] * x.N, mt.summary_stats),
                             META_AC_Allele2=all_and_leave_one_out(mt.summary_stats.AC_Allele2, mt.pheno_data.pop),
                             META_N=all_and_leave_one_out(mt.summary_stats.N, mt.pheno_data.pop))
    mt = mt.annotate_entries(META_N_pops=all_and_leave_one_out(mt.variant_exists, mt.pheno_data.pop),
                             META_AF_Allele2=mt.META_AC_Allele2 / mt.META_N,
                             META_AF_Cases=all_and_leave_one_out(mt.af_cases, mt.pheno_data.pop),
                             META_AF_Controls=all_and_leave_one_out(mt.af_controls, mt.pheno_data.pop))
    mt = mt.drop('variant_exists', 'af_cases', 'af_controls')

    # Format everything into array<struct>
    mt = mt.transmute_entries(meta_analysis=hl.map(
        lambda x: hl.struct(BETA=x[0],
                            SE=x[1],
                            Pvalue=x[2],
                            N=x[3],
                            N_pops=x[4],
                            AC_Allele2=x[5],
                            AF_Allele2=x[6],
                            AF_Cases=x[7],
                            AF_Controls=x[8]),
        hl.zip(mt.META_BETA, mt.META_SE, mt.META_Pvalue, mt.META_N, mt.META_N_pops, mt.META_AC_Allele2,
               mt.META_AF_Allele2, mt.META_AF_Cases, mt.META_AF_Controls)))

    mt = mt.annotate_cols(meta_analysis_data=hl.map(
        lambda x: hl.struct(
            n_cases=x[0],
            n_controls=x[1],
            n_cases_both_sexes=x[2],
            n_cases_females=x[3],
            n_cases_males=x[4],
            # data_type=mt.pheno_data.data_type[0],
            # meaning=mt.pheno_data.meaning[0],
            # path=mt.pheno_data.path[0],
            pop=x[5]),
        hl.zip(
            all_and_leave_one_out(mt.pheno_data.n_cases, mt.pheno_data.pop),
            all_and_leave_one_out(mt.pheno_data.n_controls, mt.pheno_data.pop),
            all_and_leave_one_out(mt.pheno_data.n_cases_both_sexes, mt.pheno_data.pop),
            all_and_leave_one_out(mt.pheno_data.n_cases_females, mt.pheno_data.pop),
            all_and_leave_one_out(mt.pheno_data.n_cases_males, mt.pheno_data.pop),
            all_and_leave_one_out(
                mt.pheno_data.pop,
                mt.pheno_data.pop,
                all_f=lambda x: x,
                loo_f=lambda i, x: hl.filter(lambda y: y != x[i], x),
            ))))

    mt.describe()
    mt.write('gs://ukb-diverse-pops/combined_results/meta_analysis.mt', overwrite=args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--overwrite', help='Overwrite data', action='store_true')
    args = parser.parse_args()

    main(args)
