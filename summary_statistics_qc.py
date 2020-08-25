#!/usr/bin/env python3

__author__ = 'konradk'

import argparse
from pprint import pprint

import gnomad.resources.grch37.gnomad as gnomad
from gnomad.utils import slack
from ukbb_pan_ancestry import *


def generate_qc_lambdas(overwrite):
    # This function could eventually read the full MT, but runs out of memory even for just EUR, so leaving split for now
    for pop in POPS:
        if pop != 'EUR':
            mt = hl.read_matrix_table(get_variant_results_path(pop, 'mt'))
            ht = generate_lambda_ht_by_freq(mt)
        else:
            # _intervals = get_n_even_intervals(10000)
            mt = hl.read_matrix_table(get_variant_results_path(pop, 'mt')).add_col_index()
            n_sections = 4
            mt = mt.annotate_cols(section=hl.find(lambda x: mt.col_idx / mt.count_cols() >= x / n_sections,
                                                  hl.range(n_sections - 1, -1, -1)))
            hts = []
            for section in range(n_sections):
                mt_section = mt.filter_cols(mt.section == section).drop('col_idx', 'section')
                ht = generate_lambda_ht_by_freq(mt_section)
                hts.append(ht.checkpoint(get_analysis_data_path('lambda', f'lambdas_part{section}', pop, 'ht'), overwrite=overwrite, _read_if_exists=not overwrite))
            ht = hts[0].union(*hts[1:])
        ht = ht.checkpoint(get_analysis_data_path('lambda', 'lambdas', pop, 'ht'), overwrite=overwrite, _read_if_exists=not overwrite)
        explode_lambda_ht(ht).export(get_analysis_data_path('lambda', 'lambdas_by_ac', pop))
        explode_lambda_ht(ht, 'af').export(get_analysis_data_path('lambda', 'lambdas_by_af', pop))


def generate_final_lambdas(overwrite):
    mt = hl.read_matrix_table(get_variant_results_path('full', 'mt'))
    qual_ht = hl.read_table(get_variant_results_qc_path())
    mt = mt.annotate_rows(**qual_ht[mt.row_key])
    ht = mt.annotate_cols(
        pheno_data=hl.zip(mt.pheno_data, hl.agg.array_agg(
            lambda ss: hl.agg.filter(~ss.low_confidence & mt.high_quality,
                hl.struct(lambda_gc=hl.methods.statgen._lambda_gc_agg(ss.Pvalue),
                          n_variants=hl.agg.count_where(hl.is_defined(ss.Pvalue)),
                          n_sig_variants=hl.agg.count_where(ss.Pvalue < 5e-8))),
            mt.summary_stats)).map(lambda x: x[0].annotate(**x[1]))
    ).cols()
    ht = ht.checkpoint(get_analysis_data_path('lambda', 'lambdas', 'full', 'ht'), overwrite=overwrite, _read_if_exists=not overwrite)
    ht.explode('pheno_data').flatten().export(get_analysis_data_path('lambda', 'lambdas', 'full', 'txt.bgz'))


def prepare_gnomad_stats(overwrite):
    ht = hl.read_table(get_ukb_vep_path())
    gwased_variants_ht = hl.read_matrix_table(get_variant_results_path('full')).rows()
    ht = ht.filter(hl.is_defined(gwased_variants_ht[ht.key]))
    datasets = ('exomes', 'genomes')
    gnomad_hts = {dataset: gnomad.public_release(dataset).ht() for dataset in datasets}
    gnomad_index_dicts = {k: gnomad_ht.index_globals().freq_index_dict for k, gnomad_ht in gnomad_hts.items()}
    gnomads = {k: gnomad_ht[ht.key] for k, gnomad_ht in gnomad_hts.items()}
    ht = ht.annotate(
        freq=[hl.struct(
            pop=pop, ac=ht.af[pop] * ht.an[pop], af=ht.af[pop] / 2, an=ht.an[pop] * 2,
            **{f'gnomad_{dataset}_{metric.lower()}':
                   gnomads[dataset].freq[gnomad_index_dicts[dataset].get(f'gnomad_{UKB_GNOMAD_POP_MAPPING[pop]}')][
                       metric]
               for dataset in datasets for metric in ('AC', 'AF', 'AN')}
        ) for pop in POPS],
        **{f'pass_gnomad_{dataset}': hl.len(gnomads[dataset].filters) == 0 for dataset in datasets}
    )
    ht.naive_coalesce(1000).write(get_analysis_data_path('gnomad_comparison', 'qc', 'full', 'ht'), overwrite=overwrite)


def generate_gnomad_plot_data():
    ht_all = hl.read_table(get_analysis_data_path('gnomad_comparison', 'qc', 'full', 'ht'))
    calling_regions = hl.import_locus_intervals('gs://gcp-public-data--broad-references/hg19/v0/wgs_calling_regions.v1.interval_list', reference_genome='GRCh37')
    genotyped_variants = hl.read_matrix_table('gs://ukb31063/ukb31063.genotype.mt').rows()
    hwe_ht = hl.import_table('gs://ukb-diverse-pops-public/misc/konradk/hwe_ldvars.tsv.gz', force_bgz=True)
    hwe_ht = hwe_ht.key_by(locus=hl.parse_locus(hwe_ht.locus, 'GRCh37'),
                           alleles=hl.parse_json(hwe_ht.alleles, hl.tarray(hl.tstr)))

    ht_all = ht_all.annotate(pass_status=hl.case(missing_false=True)
                             .when(ht_all.pass_gnomad_genomes & ht_all.pass_gnomad_exomes, 'PASS in genomes and exomes')
                             .when(ht_all.pass_gnomad_genomes | ht_all.pass_gnomad_exomes, 'PASS in at least one')
                             .when(~ht_all.pass_gnomad_genomes | ~ht_all.pass_gnomad_exomes, 'Present, never PASS')
                             .default('Not in gnomAD'),
                             in_gnomad_genomes=hl.is_defined(ht_all.pass_gnomad_genomes),
                             in_gnomad_calling=hl.str(hl.is_defined(calling_regions[ht_all.locus])),
                             in_ld_with_hwe_problematic_variant=hl.is_defined(hwe_ht[ht_all.key]),
                             genotyped_variant=hl.is_defined(genotyped_variants[ht_all.key]))

    # Variants found in gnomAD genomes have ti/tv of ~2.1, ones not in gnomAD ~1.6
    titv = ht_all.group_by('in_gnomad_genomes').aggregate(
        tis=hl.agg.count_where(hl.is_transition(ht_all.alleles[0], ht_all.alleles[1])),
        tvs=hl.agg.count_where(hl.is_transversion(ht_all.alleles[0], ht_all.alleles[1]))
    )
    titv.annotate(titv=titv.tis/titv.tvs).show()

    fc_cutoff = 2
    ht_grouped = ht_all.group_by(chrom=ht_all.locus.contig, kb=hl.int(ht_all.locus.position / 1000),
                                 gnomad_pass=ht_all.pass_gnomad_genomes).aggregate(
        res=hl.agg.array_agg(lambda x: (hl.agg.take(x.pop, 1)[0],
                                        hl.agg.count_where(x.ac >= 20),
                                        hl.agg.count_where((x.ac >= 20) & (x.af / x.gnomad_genomes_af > fc_cutoff))),
                             ht_all.freq)
    )
    ht_grouped = ht_grouped.explode('res')
    ht_grouped.select(pop=ht_grouped.res[0], n_variants=ht_grouped.res[1], n_higher_freq=ht_grouped.res[2]).export(
        get_analysis_data_path('gnomad_comparison', 'n_variants_by_region_gnomad_status', 'full'))

    for i, pop in enumerate(POPS):
        if pop in ('CSA', 'MID'): continue
        ht = ht_all.filter(ht_all.freq[i].ac >= 20)
        freq = ht.freq[i]
        ht.group_by(
            'genotyped_variant', 'pass_gnomad_genomes', 'in_gnomad_calling', 'in_ld_with_hwe_problematic_variant',
            ac0=freq.gnomad_genomes_ac == 0, freq_bin=hl.int(hl.int(hl.log10(freq.af) * 10) / 10)
        ).aggregate(n_variants=hl.agg.count()).export(get_analysis_data_path('gnomad_comparison', 'data_summary', pop))

        # def get_bounds(x):
        #     return hl.agg.filter(x > 0, hl.agg.stats(-hl.log10(x)))
        # pprint(ht.aggregate(hl.struct(gnomad=get_bounds(freq.gnomad_genomes_af), ukb=get_bounds(freq.af))))

        ht2 = hl.plot.plots._generate_hist2d_data(
            x=hl.log10(freq.gnomad_genomes_af),
            y=hl.log10(freq.af),
            range=((-4, 0), (-8, 0)), bins=(50, 100)).key_by()
        ht2.select(
            gnomad_freq=ht2.x, ukb_freq=ht2.y, n_variants=ht2.c
        ).export(get_analysis_data_path('gnomad_comparison', 'gnomad_af_v_ukb_af_hist2d', pop))

        ht2 = downsample_table_by_x_y(ht, hl.log10(freq.gnomad_genomes_af), hl.log10(freq.af),
                                      {'pass_gnomad_genomes': hl.str(ht.pass_gnomad_genomes),
                                       'genotyped_variant': hl.str(ht.genotyped_variant),
                                       'in_ld_with_hwe_problematic_variant': hl.str(ht.in_ld_with_hwe_problematic_variant),
                                       'chrom': ht.locus.contig, 'pos': hl.str(ht.locus.position)},
                                      x_field_name='gnomad_freq', y_field_name='ukb_freq')
        ht2.export(get_analysis_data_path('gnomad_comparison', 'gnomad_af_v_ukb_af_scatter', pop))

        fold_change = hl.log10(freq.af / freq.gnomad_genomes_af)
        cutoff = 3
        fold_change_filt = (hl.case()
                            .when(hl.is_nan(fold_change), hl.null(hl.tfloat64))
                            .when(fold_change >= cutoff, cutoff)
                            .when(fold_change <= -cutoff, -cutoff)
                            .default(fold_change))
        pprint(ht.aggregate(hl.agg.counter(hl.is_defined(fold_change_filt))))
        ht2 = approx_cdf_as_table(ht, fold_change_filt)
        ht2.export(get_analysis_data_path('gnomad_comparison', 'fold_change_ukb_over_gnomad_af', pop))

        n_bins = 100
        factor = n_bins / cutoff / 2
        ht_grouped = ht.group_by(fold_change_bin=hl.int(fold_change_filt * factor),
                                 freq_bin=hl.int(hl.int(hl.log10(freq.af) * 10) / 10)).aggregate(
            tis=hl.agg.count_where(hl.is_transition(ht.alleles[0], ht.alleles[1])),
            tvs=hl.agg.count_where(hl.is_transversion(ht.alleles[0], ht.alleles[1]))
        )
        ht_grouped = ht_grouped.key_by()
        ht_grouped = ht_grouped.annotate(fold_change_bin=ht_grouped.fold_change_bin / factor)
        ht_grouped.export(get_analysis_data_path('gnomad_comparison', 'titv_by_fold_change_ukb_over_gnomad_af', pop))


        ctt = hl.chi_squared_test(hl.int(freq.ac), hl.int(freq.an - freq.ac), freq.gnomad_genomes_ac,
                                  freq.gnomad_genomes_an - freq.gnomad_genomes_ac)
        # hl.plot.show(hl.plot.scatter(hl.log10(ctt.odds_ratio), -hl.log10(ctt.p_value),
        #                              {'pass_gnomad_genomes': hl.str(ht.pass_gnomad_genomes)
        #                               }))

        ht2 = hl.plot.plots._generate_hist2d_data(
            x=hl.if_else(hl.log10(ctt.odds_ratio) > 3, 3, hl.log10(ctt.odds_ratio)),
            y=hl.if_else(-hl.log10(ctt.p_value) > 10, 10, -hl.log10(ctt.p_value)),
            range=((-3, 3), (0, 10)), bins=100).key_by()
        ht2.select(
            odds_ratio=ht2.x, p_value=ht2.y, n_variants=ht2.c
        ).export(get_analysis_data_path('gnomad_comparison', 'or_vs_p_hist2d', pop))

    chisq_results = ht_all.freq.filter(lambda x: filter_by_chisq(x)).map(lambda x: x.pop)
    ht_grouped = ht_all.group_by('in_gnomad_genomes', failing_pops=hl.delimit(chisq_results)).aggregate(
        n=hl.agg.count(),
        res=hl.agg.array_agg(lambda x: (hl.agg.take(x.pop, 1)[0], hl.agg.count_where(x.ac >= 20)),
                             ht_all.freq)
    )
    ht_grouped.select('n').export(get_analysis_data_path('gnomad_comparison', 'n_variants_by_failing_pops', 'full'))
    ht_grouped = ht_grouped.explode('res')
    ht_grouped.select(pop=ht_grouped.res[0], n_variants=ht_grouped.res[1]).export(
        get_analysis_data_path('gnomad_comparison', 'n_variants_by_failing_pops_by_pop', 'full'))


def approx_cdf_as_table(ht, field):
    res = ht.aggregate(hl.agg.approx_cdf(field), _localize=False)
    ht2 = hl.utils.range_table(1).annotate(data=hl.zip(res.values, res.ranks[1:])).explode('data')
    return ht2.select(values=ht2.data[0], ranks=ht2.data[1]).key_by().drop('idx')


def filter_by_chisq(freq, or_cutoff: float = 2.0, p_cutoff: float = 1e-6):
    csq = hl.chi_squared_test(hl.int(freq.ac), hl.int(freq.an - freq.ac), freq.gnomad_genomes_ac,
                        freq.gnomad_genomes_an - freq.gnomad_genomes_ac)
    return (csq.odds_ratio > or_cutoff) & (csq.p_value < p_cutoff)


def generate_variant_qc_file(overwrite: bool = False):
    ht = hl.read_table(get_analysis_data_path('gnomad_comparison', 'qc', 'full', 'ht'))
    chisq_results = ht.freq.filter(lambda x: hl.is_defined(x.gnomad_genomes_af) & ~filter_by_chisq(x))
    ht = ht.annotate(
        n_passing_populations=hl.len(chisq_results),
        high_quality=ht.pass_gnomad_genomes & (hl.len(chisq_results) == 4)
    ).drop('af', 'an')
    ht = annotate_nearest_gene(ht, add_only_gene_symbols_as_str=True)
    info_score = hl.read_table(ukb_imputed_info_ht_path)
    ht = ht.annotate(info=info_score[ht.key].info).checkpoint(get_variant_results_qc_path(), overwrite, _read_if_exists=not overwrite)

    ht = ht.drop('vep', 'pass_gnomad_exomes')
    def get_freq_data(pop):
        metrics = ('ac', 'an', 'af')
        freq = ht.freq.find(lambda x: x.pop == pop).drop('pop', *[f'gnomad_exomes_{m}' for m in metrics])
        if pop in ('CSA', 'MID'):
            freq = freq.drop(*[f'gnomad_genomes_{m}' for m in metrics])
        return freq
    ht2 = ht.transmute(**{pop: get_freq_data(pop) for pop in POPS}).flatten()
    def rename_pop_metric(pop_metric):
        if '.' not in pop_metric: return pop_metric
        pop, metric = pop_metric.split('.')
        return f'{metric}_{pop}'
    ht2 = ht2.rename({x: rename_pop_metric(x) for x in ht2.row_value})
    ht = locus_alleles_to_chr_pos_ref_alt(ht2, True)
    ht.export(get_variant_results_qc_path('txt.bgz'))


def main(args):
    hl.init(default_reference='GRCh37', log='/combine_results.log', branching_factor=8)

    if args.compute_qc_lambdas:
        generate_qc_lambdas(args.overwrite)

    if args.compute_variant_metrics:
        prepare_gnomad_stats(args.overwrite)
        generate_gnomad_plot_data()
        generate_variant_qc_file(args.overwrite)

    if args.compute_final_lambdas:
        generate_final_lambdas(args.overwrite)

    if args.gwas_run_stats:
        mt = hl.read_matrix_table(get_variant_results_path('full'))
        ht = mt.cols()
        pprint(f'Phenos by number of pops: {ht.aggregate(hl.agg.counter(hl.len(ht.pheno_data)))}')
        ht = ht.explode('pheno_data')
        pprint(f'Total GWAS: {ht.count()}')

        ht.group_by('trait_type', pop=ht.pheno_data.pop).aggregate(n_phenos=hl.agg.count()).show()
        trait_types = ht.aggregate(hl.agg.collect_as_set(ht.trait_type))
        ht.group_by(pop=ht.pheno_data.pop).aggregate(
            n_samples=hl.agg.max(ht.pheno_data.n_cases + hl.or_else(ht.pheno_data.n_controls, 0)),
            total_n_phenos=hl.agg.count(),
            **{trait_type: hl.agg.count_where(ht.trait_type == trait_type) for trait_type in trait_types}).show()

        # Confirmed that the overall number matches those with genotype+phenotype data:
        # genotyped_samples = get_filtered_mt('22').cols()
        # ht = hl.read_matrix_table(get_ukb_pheno_mt_path()).rows()
        # ht = ht.key_by(s=hl.str(ht.userId))
        # all_samples = genotyped_samples.filter(hl.is_defined(ht[genotyped_samples.key]))
        # all_samples.aggregate(hl.agg.counter(all_samples.pop))

        # x = ht.filter(ht.trait_type == 'biomarkers')
        # x.group_by(x.phenocode, x.description).aggregate(n=hl.agg.count()).show(40)

        # Total number of statistical tests
        lambda_ht = hl.read_table(get_analysis_data_path('lambda', 'lambdas', 'full', 'ht')).explode('pheno_data')
        lambda_ht.aggregate(hl.agg.sum(lambda_ht.pheno_data.n_variants))
        lambda_ht.group_by(lambda_ht.pheno_data.pop).aggregate(total_tests=hl.agg.sum(lambda_ht.pheno_data.n_variants)).show()

    if args.release_stats_qc:
        ht = hl.read_table(get_analysis_data_path('lambda', 'lambdas', 'full', 'ht')).explode('pheno_data')
        ht = ht.filter(filter_lambda_gc(ht.pheno_data.lambda_gc))

        # Total number of significant hits
        res = ht.aggregate(hl.struct(total_sig_variants=hl.agg.sum(ht.pheno_data.n_sig_variants),
                                     total_variants=hl.agg.sum(ht.pheno_data.n_variants)))
        print(f'Got {res.total_sig_variants:,} significant hits (out of {res.total_variants:,} total tests)')



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--compute_variant_metrics', help='Overwrite everything', action='store_true')
    parser.add_argument('--compute_qc_lambdas', help='Overwrite everything', action='store_true')
    parser.add_argument('--compute_final_lambdas', help='Overwrite everything', action='store_true')
    parser.add_argument('--gwas_run_stats', help='Overwrite everything', action='store_true')
    parser.add_argument('--release_stats_qc', help='Overwrite everything', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        from slack_token_pkg.slack_creds import slack_token
        with slack.slack_notifications(slack_token, args.slack_channel):
            main(args)
