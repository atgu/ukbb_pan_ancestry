#!/usr/bin/env python3

__author__ = 'konradk'

import sys
import argparse
import logging
logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s", level='INFO', filename='saige_pipeline.log')

from ukbb_common import *
from ukbb_pan_ancestry import *

logger = logging.getLogger("saige_pan_ancestry")
logger.addHandler(logging.StreamHandler(sys.stderr))


def main(args):
    hl.init(log=f'/tmp/saige_temp_hail.log')
    hts = []
    for pop in POPS:
        ht = hl.import_table(get_results_timing_tsv_path('saige', pop), types={'saige_time': hl.tfloat64}, min_partitions=100)
        ht = ht.annotate(locus=hl.parse_locus(ht.locus))
        ht = ht.key_by('pop', 'trait_type', 'pheno', 'coding', 'locus')
        hts.append(ht)
    ht = hts[0].union(*hts[1:])
    ht.write(get_results_timing_ht_path('saige'), True)

    ht = hl.import_table(get_results_timing_tsv_path('null_model'), impute=True)
    ht = ht.key_by('pop', trait_type=ht.trait, pheno=ht.pheno, coding=ht.coding)
    ht.write(get_results_timing_ht_path('null_model'), True)

    ht = hl.read_table(get_results_timing_ht_path('saige'))
    ht = ht.group_by('pop', 'trait_type', 'pheno', 'coding').aggregate(saige_time=hl.agg.sum(ht.saige_time))
    null_ht = hl.read_table(get_results_timing_ht_path('null_model'))
    ht.annotate(**null_ht[ht.key]).write(get_results_timing_ht_path('full'), True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
