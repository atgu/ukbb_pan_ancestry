#!/usr/bin/env python3

__author__ = 'konradk'

import sys
import argparse
import logging
logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s", level='INFO', filename='saige_pipeline.log')

from ukb_common import *
import time

from ukbb_pan_ancestry import *
from ukb_common.utils.saige_pipeline import *

logger = logging.getLogger("saige_pan_ancestry")
logger.addHandler(logging.StreamHandler(sys.stderr))
bucket = 'gs://ukb-diverse-pops'
root = f'{bucket}/results'

HAIL_DOCKER_IMAGE = 'gcr.io/ukbb-diversepops-neale/hail_utils:3.5'
SAIGE_DOCKER_IMAGE = 'wzhou88/saige:0.36.3'
QQ_DOCKER_IMAGE = 'konradjk/saige_qq:0.2'


def get_phenos_to_run(pop: str):
    ht = hl.read_table(get_phenotype_summary_path('full'))
    ht = ht.filter(ht.pop == pop)
    min_cases = MIN_CASES_EUR if pop == 'EUR' else MIN_CASES
    ht = ht.filter(~hl.set({'raw', 'icd9'}).contains(ht.coding) &
                   ~hl.set({'22601', '22617', '20024', '41230', '41210'}).contains(ht.pheno) &
                   (ht.n_cases_by_pop >= min_cases) &
                   (ht.n_cases_both_sexes >= MIN_CASES_ALL))

    fields = ('pheno', 'coding', 'data_type')
    output = set([tuple(x[field] for field in fields) for x in ht.key_by().select(*fields).collect()])
    output = output.union(PILOT_PHENOTYPES)
    return output


def format_pheno_dir(pheno):
    return pheno.replace("/", "_")


def get_inverse_normalize_status(null_glmm_log):
    status = 'Unknown'
    with hl.hadoop_open(null_glmm_log) as f:
        for line in f:
            if line.startswith('$invNormalize'):
                try:
                    status = f.readline().strip().split()[1]
                except:
                    logger.warn(f'Could not load inv_norm status from {line} in {null_glmm_log}.')
    return status.capitalize()


def main(args):
    hl.init(log='/tmp/saige_temp_hail.log')

    output = 'pop\tpheno\tcoding\ttrait\theritability\tinv_normalized\n'
    for pop in POPS:
        logger.info(f'Setting up {pop}...')
        phenos_to_run = get_phenos_to_run(pop)

        null_model_dir = f'{root}/null_glmm/{pop}'
        for i, pheno_coding_trait in enumerate(phenos_to_run):
            pheno, coding, trait_type = pheno_coding_trait
            quantitative_trait = trait_type in ('continuous', 'biomarkers')
            null_glmm_log = f'{null_model_dir}/{trait_type}-{format_pheno_dir(pheno)}-{coding}.variant.log'
            pct = "\t".join(pheno_coding_trait)
            try:
                heritability = get_heritability_from_log(null_glmm_log, quantitative_trait)
                inv_normalized = get_inverse_normalize_status(null_glmm_log)
                output += f'{pop}\t{pct}\t{heritability}\t{inv_normalized}\n'
            except:
                logger.info(f'Missing: {null_glmm_log}')

    with hl.hadoop_open(get_heritability_txt_path(), 'w') as f:
        f.write(output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite_pheno_data', help='Run single variant SAIGE', action='store_true')
    args = parser.parse_args()

    try_slack('@konradjk', main, args)

