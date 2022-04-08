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
bucket = 'gs://ukb-diverse-pops'
root = f'{bucket}/results'


def format_pheno_dir(pheno):
    return pheno.replace("/", "_")


def main(args):
    hl.init(log='/tmp/saige_temp_hail.log')
    log_extension = '.variant.log'

    output = 'pop\ttrait\tpheno\tcoding\tnull_model_time\tnull_model_wall\n'
    for pop in POPS:
        logger.info(f'Setting up {pop}...')
        null_model_dir = f'{root}/null_glmm/{pop}'
        all_phenos_dir = hl.hadoop_ls(null_model_dir)
        all_phenos_dir = [x['path'] for x in all_phenos_dir if x['path'].endswith(log_extension)]
        logger.info(f'Got {len(all_phenos_dir)} phenotypes for {pop}...')
        for fname in all_phenos_dir:
            try:
                trait_pheno_coding = fname.split(log_extension)[0].split('/')[-1]
                trait_type, pheno_coding = trait_pheno_coding.split('-', 1)
                pheno, coding = pheno_coding.rsplit('-', 1)
            except:
                logger.warning(f'Malformed filename at {fname}')
                continue

            null_model_time, null_model_wall = get_null_model_timing(fname)

            tpc = "\t".join([trait_type, pheno, coding])
            logger.info(f'{pop}\t{tpc}\t{null_model_time}\t{null_model_wall}')

            output += f'{pop}\t{tpc}\t{null_model_time}\t{null_model_wall}\n'

    with hl.hadoop_open(f'{root}/misc/timings_null_model.txt', 'w') as f:
        f.write(output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite_pheno_data', help='Run single variant SAIGE', action='store_true')
    args = parser.parse_args()

    try_slack('@konradjk', main, args)

