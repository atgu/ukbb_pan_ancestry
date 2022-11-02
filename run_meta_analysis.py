#!/usr/bin/env python
# coding: utf-8
import argparse
import hail as hl
from ukbb_pan_ancestry.resources.results import get_meta_analysis_results_path
from ukbb_pan_ancestry.utils.results import load_final_sumstats_mt
from ukbb_pan_ancestry.utils.meta_analysis import run_meta_analysis


def main(args):
    # Read in all sumstats
    mt = load_final_sumstats_mt(
        filter_phenos=True,
        filter_variants=False,
        filter_sumstats=(not args.keep_low_confidence_variants),
        separate_columns_by_pop=False,
        annotate_with_nearest_gene=False,
        filter_pheno_h2_qc=(not args.keep_low_heritability_pairs),
    )
    mt = run_meta_analysis(mt)

    mt.describe()
    mt.write(
        get_meta_analysis_results_path(filter_pheno_h2_qc=(not args.keep_low_heritability_pairs)),
        overwrite=args.overwrite,
    )

    hl.copy_log("gs://ukb-diverse-pops/combined_results/meta_analysis.log")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--keep-low-confidence-variants",
        help="Keep variants with low confidence flag for meta-analysis",
        action="store_true",
    )
    parser.add_argument(
        "--keep-low-heritability-pairs",
        help="Keep trait-ancestry pairs that failed at heritability QC",
        action="store_true",
    )
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    args = parser.parse_args()

    main(args)
