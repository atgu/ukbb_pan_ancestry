#!/usr/bin/env python3

__author__ = 'konradk'

import argparse
from pprint import pprint

from gnomad.utils import slack
from ukb_common import *
from ukbb_pan_ancestry import *

def remove_phenos_from_analysis(mt: hl.MatrixTable):
    pass


def annotate_mt_with_largest_meta_analysis(mt):
    meta_mt = hl.read_matrix_table(get_meta_analysis_results_path())
    meta_analysis_join = meta_mt[mt.row_key, mt.col_key].meta_analysis
    return mt.annotate_entries(meta=hl.or_missing(hl.len(meta_analysis_join) > 0, meta_analysis_join[0]))


def main(args):
    hl.init(default_reference='GRCh37')

    if args.write_gene_intervals:
        create_genome_intervals_file().write(get_gene_intervals_path('GRCh37'), args.overwrite)

    if args.compute_top_p:
        mt = load_final_sumstats_mt(separate_columns_by_pop=False, add_only_gene_symbols_as_str=True)
        mt = annotate_mt_with_largest_meta_analysis(mt)

        sumstats_with_pop = hl.zip_with_index(mt.summary_stats).map(
            lambda x: x[1].annotate(pop=mt.pheno_data[x[0]].pop)
        )
        sumstats_with_pop_pheno_meta = sumstats_with_pop.map(
            lambda x: x.annotate(meta=mt.meta, **mt.col_key, **{x: mt[x] for x in PHENO_DESCRIPTION_FIELDS})
        )
        sumstats_with_pheno_meta_by_pop = [sumstats_with_pop_pheno_meta.find(lambda ss: ss.pop == pop) for pop in POPS]

        meta_with_pheno_and_per_pop = mt.meta.annotate(
            pop_data=sumstats_with_pop,
            **mt.col_key, **{x: mt[x] for x in PHENO_DESCRIPTION_FIELDS})

        def agg_top_p(ss_pop, significant_only: bool = False):
            p = -hl.abs(ss_pop.BETA / ss_pop.SE)  # ss_pop.Pvalue
            top_ss = hl.agg.take(ss_pop, 1, ordering=p)
            p_filter = hl.is_defined(p) & ~hl.is_nan(p)
            if significant_only:
                p_filter &= (p < 5e-8)
            return hl.agg.filter(p_filter, hl.or_missing(hl.len(top_ss) > 0, top_ss[0]))

        ht = mt.annotate_rows(
            top_p=[agg_top_p(ss_pop) for ss_pop in sumstats_with_pheno_meta_by_pop],
            top_meta_p=agg_top_p(meta_with_pheno_and_per_pop)
        ).rows().naive_coalesce(1000).drop('vep', 'freq')
        ht.write(get_analysis_data_path('sig_hits', 'top_p_by_variant', 'full', 'ht'), overwrite=args.overwrite)

    if args.export_top_p:
        ht_full = hl.read_table(get_analysis_data_path('sig_hits', 'top_p_by_variant', 'full', 'ht'))

        ht = locus_alleles_to_chr_pos_ref_alt(ht_full.annotate(global_position=ht_full.locus.global_position()), True)
        ht = ht.transmute(**ht.top_meta_p).drop('top_p')
        ht = ht.transmute(eur_data=ht.pop_data.find(lambda x: x.pop == 'EUR'))

        # Downsample in EUR p vs meta p partial-log-log space
        def partial_log_log(p):
            neglog_p = -hl.log10(p)
            return hl.if_else(neglog_p > 10, 10 * hl.log10(neglog_p), neglog_p)

        # ht2 = downsample_table_by_x_y(ht, partial_log_log(ht.eur_data.Pvalue), partial_log_log(ht.Pvalue),
        #                               {'nearest_genes': ht.nearest_genes,
        #                                **{x: hl.str(ht[x]) for x in ('chrom', 'pos', 'ref', 'alt') + PHENO_KEY_FIELDS + PHENO_DESCRIPTION_FIELDS}},
        #                               x_field_name='eur_pvalue_ll', y_field_name='meta_pvalue_ll', n_divisions=500)
        # ht2.export(get_analysis_data_path('sig_hits', 'eur_p_vs_meta_p', 'EUR'))

        ht2 = downsample_table_by_x_y(ht, ht.global_position, partial_log_log(ht.Pvalue),
                                      {'nearest_genes': ht.nearest_genes,
                                       **{x: hl.str(ht[x]) for x in ('chrom', 'pos', 'ref', 'alt') + PHENO_KEY_FIELDS + PHENO_DESCRIPTION_FIELDS}},
                                      x_field_name='global_position', y_field_name='meta_pvalue_ll', n_divisions=500)
        ht2.export(get_analysis_data_path('sig_hits', 'downsampled_manhattan', 'meta'))
        return
        # Get variants that are 5e-8 in meta, but not in EUR
        ht = ht.filter((ht.eur_data.Pvalue > 5e-8) & (ht.Pvalue < 5e-8)).flatten()
        ht.export(get_analysis_data_path('sig_hits', 'newly_significant_variants', 'full'))

        ht_full = locus_alleles_to_chr_pos_ref_alt(ht_full.filter(ht_full.locus.contig == '2'), True)
        # Top meta p value by variant
        ht = ht_full.filter(ht_full.top_meta_p.Pvalue < 0.01)
        ht.transmute(**ht.top_meta_p).drop('top_p', 'pop_data').export(
            get_analysis_data_path('sig_hits', 'top_meta_p_by_variant', 'full'))

        # Top meta p value by variant with data from each population
        ht.transmute(**ht.top_meta_p).drop('top_p').explode('pop_data').flatten().export(
            get_analysis_data_path('sig_hits', 'top_meta_p_by_variant_with_pop', 'full'))

        # Top p value for each population (with meta)
        ht = ht_full.explode('top_p')
        ht = ht.filter(ht.top_p.Pvalue < 0.01)
        ht.transmute(**ht.top_p.annotate(meta_for_top_pheno_for_pop=ht.top_p.meta).drop('meta', 'description_more'),
                     **{f'meta_{k}': v for k, v in ht.top_meta_p.items() if k != 'description_more'}).flatten().export(
            get_analysis_data_path('sig_hits', 'top_p_by_variant_by_pop', 'full'))

    if args.compute_sig_pops:
        mt = load_final_sumstats_mt(separate_columns_by_pop=False)
        mt = mt.annotate_entries(
            sig_pops=get_sig_pops(mt)
        )
        mt = mt.annotate_cols(
            sig_pops_by_pheno=hl.agg.counter(hl.delimit(hl.sorted(mt.sig_pops)))
        )
        ht = mt.cols().checkpoint(get_analysis_data_path('sig_hits', 'sig_hits_pops_by_pheno', 'full', 'ht'), overwrite=True)
        ht = ht.transmute(output=hl.zip(ht.sig_pops_by_pheno.keys(), ht.sig_pops_by_pheno.values())).explode('output')
        ht.transmute(pop_group=ht.output[0], count=ht.output[1]).drop('pheno_data', 'pheno_indices').export(get_analysis_data_path('sig_hits', 'sig_hits_pops_by_pheno', 'full'))

    if args.hits_by_pheno:
        mt = load_final_sumstats_mt(separate_columns_by_pop=False, add_only_gene_symbols_as_str=True)
        mt = mt.annotate_entries(
            sig_pops=get_sig_pops(mt)
        )
        mt = mt.group_rows_by(
            contig=mt.locus.contig,
            genes=mt.nearest_genes,
        ).partition_hint(100).aggregate(
            sig_pops=hl.agg.explode(lambda x: hl.agg.collect_as_set(x), mt.sig_pops)
        )
        ht = mt.filter_entries(hl.len(mt.sig_pops) > 0).drop('pheno_data', 'pheno_indices').entries()
        ht = ht.checkpoint(get_analysis_data_path('sig_hits', 'sig_hits_by_pheno', 'full', 'ht'), overwrite=True)
        ht.annotate(sig_pops=hl.delimit(hl.sorted(hl.array(ht.sig_pops)))).export(
            get_analysis_data_path('sig_hits', 'sig_hits_by_pheno', 'full')
        )

    if args.meta_analysis_hits:
        mt = load_final_sumstats_mt(separate_columns_by_pop=False, add_only_gene_symbols_as_str=True)#, load_contig='22')
        mt = annotate_mt_with_largest_meta_analysis(mt)
        mt = mt.select_entries(sig_meta=hl.or_else(mt.meta.Pvalue < 5e-8, False),
                               sig_eur=hl.or_else(get_sig_pops(mt).contains('EUR'), False))
        pprint(mt.aggregate_entries(hl.struct(
            significant_associations=hl.agg.counter(mt.entry),
            significant_loci=hl.len(hl.agg.filter(mt.sig_meta, hl.agg.collect_as_set(mt.row_key))),
            significant_loci_EUR=hl.len(hl.agg.filter(mt.sig_eur, hl.agg.collect_as_set(mt.row_key))),
            significant_loci_not_EUR=hl.len(hl.agg.filter(mt.sig_meta & ~mt.sig_eur, hl.agg.collect_as_set(mt.row_key))),
            significant_traits=hl.len(hl.agg.filter(mt.sig_meta, hl.agg.collect_as_set(mt.col_key))),
            significant_traits_EUR=hl.len(hl.agg.filter(mt.sig_eur, hl.agg.collect_as_set(mt.col_key))),
            significant_traits_not_EUR = hl.len(hl.agg.filter(mt.sig_meta & ~mt.sig_eur, hl.agg.collect_as_set(mt.col_key))),
        )))
        # {'significant_associations': {Struct(sig_meta=True, sig_eur=False): 901199,
        #                               Struct(sig_meta=False, sig_eur=False): 108905474188,
        #                               Struct(sig_meta=True, sig_eur=True): 12934323,
        #                               Struct(sig_meta=False, sig_eur=True): 1056739},
        #  'significant_loci': 1790087,
        #  'significant_loci_EUR': 1826305,
        #  'significant_loci_not_EUR': 527517,
        #  'significant_traits': 3850,
        #  'significant_traits_EUR': 3752,
        #  'significant_traits_not_EUR': 1725}


def get_sig_pops(mt):
    return (hl.zip_with_index(mt.summary_stats)
            .filter(lambda x: x[1].Pvalue < 5e-8)
            .map(lambda x: mt.pheno_data[x[0]].pop))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--write_gene_intervals', help='Overwrite everything', action='store_true')
    parser.add_argument('--compute_top_p', help='Overwrite everything', action='store_true')
    parser.add_argument('--export_top_p', help='Overwrite everything', action='store_true')
    parser.add_argument('--compute_sig_pops', help='Overwrite everything', action='store_true')
    parser.add_argument('--meta_analysis_hits', help='Overwrite everything', action='store_true')
    parser.add_argument('--hits_by_pheno', help='Overwrite everything', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        from slack_token_pkg.slack_creds import slack_token
        with slack.slack_notifications(slack_token, args.slack_channel):
            main(args)
