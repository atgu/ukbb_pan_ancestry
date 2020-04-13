#!/usr/bin/env python3

__author__ = 'konradk'

import sys
import argparse
import logging
logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s", level='INFO', filename='saige_pipeline.log')

from ukb_common import *
from ukbb_pan_ancestry import *

logger = logging.getLogger("saige_pan_ancestry")
logger.addHandler(logging.StreamHandler(sys.stderr))
root = f'{bucket}/results'


def main(args):
    pop = args.pop
    hl.init(log=f'/tmp/saige_temp_hail_{pop}.log')
    log_extension = '.variant.log'
    chromosomes = list(map(str, range(1, 23))) + ['X']
    reference = 'GRCh37'
    chrom_lengths = hl.get_reference(reference).lengths

    output_path = get_results_timing_tsv_path('saige', pop)
    out_fname = output_path.rsplit('/', 1)[1]

    f = open(f'/tmp/{out_fname}', 'w')
    f.write('pop\ttrait_type\tpheno\tcoding\tlocus\tsaige_time\n')
    logger.info(f'Setting up {pop}...')

    result_dir = f'{root}/result/{pop}'
    chunk_size = int(5e6) if pop != 'EUR' else int(1e6)

    all_phenos_dir = [x['path'] for x in hl.hadoop_ls(result_dir) if x['is_dir']]
    if args.trait_type:
        all_phenos_dir = [x for x in all_phenos_dir if args.trait_type in x]
    logger.info(f'Got {len(all_phenos_dir)} phenotypes for {pop}...')

    for fname in all_phenos_dir:
        all_files = []
        trait_pheno_coding = fname.split(log_extension)[0].split('/')[-1]
        try:
            trait_type, pheno_coding = trait_pheno_coding.split('-', 1)
            pheno, coding = pheno_coding.rsplit('-', 1)
        except:
            logger.warning(f'Failed to parse: {trait_pheno_coding}')
            continue
        pheno_results_dir = f'{result_dir}/{trait_type}-{format_pheno_dir(pheno)}-{coding}'

        for chromosome in chromosomes:
            chrom_length = chrom_lengths[chromosome]
            for start_pos in range(1, chrom_length, chunk_size):
                results_path = f'{pheno_results_dir}/result_{format_pheno_dir(pheno)}_{chromosome}_{str(start_pos).zfill(9)}{log_extension}'
                all_files.append(results_path)
        logger.info(f'Got {len(all_files)} files for {pop}, pheno: {trait_pheno_coding}...')

        for locus, saige_time in get_saige_timing_grep(all_files):
            tpc = "\t".join([trait_type, pheno, coding])
            f.write(f'{pop}\t{tpc}\t{locus}\t{saige_time}\n')
    f.close()

    hl.hadoop_copy(f'file:///tmp/{out_fname}', output_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--pop')
    parser.add_argument('--trait_type')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
