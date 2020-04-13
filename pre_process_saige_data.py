#!/usr/bin/env python3

__author__ = 'konradk'

from ukb_common import *
from ukbb_pan_ancestry import *


def main(args):
    hl.init(default_reference='GRCh37')

    pops = args.pops.split(',') if args.pops else POPS

    if args.prepare_genotype_data:
        load_all_mfi_data().write(ukb_imputed_info_ht_path, args.overwrite)

    if args.genotype_summary:
        variants = hl.read_table(ukb_imputed_info_ht_path)
        print(variants.count())
        variants = variants.filter(variants.info > 0.8)
        print(variants.count())
        meta_ht = hl.import_table(get_ukb_meta_pop_tsv_path(), impute=True, types={'s': hl.tstr}, key='s')

        mt = get_ukb_imputed_data('all', variant_list=variants, entry_fields=('dosage', ))
        mt = mt.annotate_cols(**meta_ht[mt.col_key])
        ht = mt.annotate_rows(af=hl.agg.group_by(mt.pop, hl.agg.mean(mt.dosage)),
                              an=hl.agg.group_by(mt.pop, hl.agg.count_where(hl.is_defined(mt.dosage)))).rows()
        ht = ht.checkpoint(get_ukb_af_ht_path(False), args.overwrite, _read_if_exists=not args.overwrite)

        mt = get_ukb_imputed_data('X', variant_list=variants, entry_fields=('dosage', ))
        mt = mt.annotate_cols(**meta_ht[mt.col_key])
        ht_x = mt.annotate_rows(af=hl.agg.group_by(mt.pop, hl.agg.mean(mt.dosage)),
                              an=hl.agg.group_by(mt.pop, hl.agg.count_where(hl.is_defined(mt.dosage)))).rows()
        ht = ht.union(ht_x)
        ht = ht.checkpoint(get_ukb_af_ht_path(), args.overwrite, _read_if_exists=not args.overwrite)
        ht = ht.naive_coalesce(1000).checkpoint(get_ukb_af_ht_path(repart=True), args.overwrite, _read_if_exists=not args.overwrite)

        print(ht.aggregate(hl.struct(
            # hist=hl.agg.hist(hl.sum(ht.an.values()), 0, total_samples, 10),  # No missing data
            # fraction_missingness=hl.agg.fraction(hl.sum(ht.an.values()) < total_samples),
            # number_sites_above_001=hl.agg.count_where(hl.any(lambda x: x / 2 > 0.001, ht.af.values())),
            # number_sites_above_005=hl.agg.count_where(hl.any(lambda x: x / 2 > 0.005, ht.af.values())),
            # number_sites_above_01=hl.agg.count_where(hl.any(lambda x: x / 2 > 0.01, ht.af.values())),
            # number_sites_above_05=hl.agg.count_where(hl.any(lambda x: x / 2 > 0.05, ht.af.values())),
            # number_sites_above_10=hl.agg.count_where(hl.any(lambda x: x / 2 > 0.1, ht.af.values()))
            number_sites_above_mac_20=hl.agg.count_where(hl.any(lambda x: ht.af[x] * ht.an[x] >= 20, hl.literal(POPS))),
            number_run_sites_above_mac_20=hl.agg.sum(hl.sum(hl.map(lambda x: hl.int(ht.af[x] * ht.an[x] >= 20), hl.literal(POPS))))
        )))

    if args.vep:
        call_stats_ht = hl.read_table(get_ukb_af_ht_path())
        ht = vep_or_lookup_vep(call_stats_ht)
        ht.write(get_ukb_vep_path(), args.overwrite)

    if args.create_plink_file:
        window = args.window
        # Note: 1e7 LD pruning for EUR was run with r2=0.05, and chr8 inversion and HLA were removed
        r2 = 0.1 if args.window != '1e7' else 0.05
        iteration = 1
        for pop in pops:
            call_stats_ht = hl.read_table(get_ukb_af_ht_path(with_x = False))
            mt = get_filtered_mt(pop=pop, imputed=False)
            n_samples = mt.count_cols()
            print(f'Got {n_samples} samples for {pop}...')
            mt = filter_to_autosomes(mt)
            callstats = call_stats_ht[mt.row_key]
            mt = mt.filter_rows(callstats.af[pop] > 0.01)

            mt = mt.checkpoint(get_ukb_grm_mt_path(pop, iteration), _read_if_exists=not args.overwrite, overwrite=args.overwrite)
            if pop == 'EUR':
                # Common inversion taken from Table S4 of https://www.ncbi.nlm.nih.gov/pubmed/27472961
                # Also removing HLA, from https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh37
                mt = mt.filter_rows(~hl.parse_locus_interval('8:8055789-11980649').contains(mt.locus) &
                                    ~hl.parse_locus_interval('6:28477797-33448354').contains(mt.locus))
            mt = mt.unfilter_entries()
            if not args.omit_ld_prune:
                ht = hl.ld_prune(mt.GT, r2=float(r2), bp_window_size=int(float(window)))
                ht.write(get_ukb_grm_pruned_ht_path(pop, window), overwrite=args.overwrite)
            ht = hl.read_table(get_ukb_grm_pruned_ht_path(pop, window))
            mt = mt.filter_rows(hl.is_defined(ht[mt.row_key]))
            if pop == 'EUR':
                mt = mt.filter_rows(hl.rand_bool(0.55))

            if args.overwrite or not hl.hadoop_exists(f'{get_ukb_grm_plink_path(pop, iteration, window)}.bed'):
                hl.export_plink(mt, get_ukb_grm_plink_path(pop, iteration, window))

            mt = get_filtered_mt(chrom='22', pop=pop)
            if args.overwrite or not hl.hadoop_exists(get_ukb_samples_file_path(pop, iteration)):
                with hl.hadoop_open(get_ukb_samples_file_path(pop, iteration), 'w') as f:
                    f.write('\n'.join(mt.s.collect()) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--pheno_summary', action='store_true')
    parser.add_argument('--prepare_genotype_data', action='store_true')
    parser.add_argument('--genotype_summary', action='store_true')
    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--omit_ld_prune', help='Overwrite everything', action='store_true')
    parser.add_argument('--window', help='Overwrite everything', default='1e6')
    parser.add_argument('--pops', help='Comma-separated list of pops to run')
    parser.add_argument('--create_plink_file', help='Overwrite everything', action='store_true')
    parser.add_argument('--vep', help='Overwrite everything', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)


