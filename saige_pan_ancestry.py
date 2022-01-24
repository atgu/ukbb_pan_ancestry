#!/usr/bin/env python3

__author__ = 'konradk'

import sys
import copy
import argparse
import logging
logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s", level='INFO', filename='saige_pipeline.log')

from gnomad.utils import slack
from ukb_common import *
import time
import re

from ukbb_pan_ancestry import *
from ukb_common.utils.saige_pipeline import *

logger = logging.getLogger("saige_pan_ancestry")
logger.addHandler(logging.StreamHandler(sys.stderr))
bucket = 'gs://ukb-diverse-pops'
root = f'{bucket}/results'

HAIL_DOCKER_IMAGE = 'gcr.io/ukbb-diversepops-neale/hail_utils:6.1'
SAIGE_DOCKER_IMAGE = 'gcr.io/ukbb-diversepops-neale/saige:0.5'
QQ_DOCKER_IMAGE = 'konradjk/saige_qq:0.2'


def main(args):
    hl.init(log='/tmp/saige_temp_hail.log')

    # num_pcs = 20
    num_pcs = 10
    start_time = time.time()
    basic_covars = ['sex', 'age', 'age2', 'age_sex', 'age2_sex']
    covariates = ','.join(basic_covars + [f'PC{x}' for x in range(1, num_pcs + 1)])
    n_threads = 8
    analysis_type = "variant"
    chromosomes = list(map(str, range(1, 23))) + ['X']
    reference = 'GRCh37'
    chrom_lengths = hl.get_reference(reference).lengths
    iteration = 1
    pops = args.pops.split(',') if args.pops else POPS

    # if args.local_test:
    #     backend = hb.LocalBackend(gsa_key_file='/Users/konradk/.hail/ukb-diverse-pops.json')
    # else:

    backend = hb.ServiceBackend(billing_project='ukb_diverse_pops',
                                bucket=temp_bucket.split('gs://', 1)[-1])
    for pop in pops:
        p = hb.Batch(name=f'saige_pan_ancestry_{pop}', backend=backend, default_image=SAIGE_DOCKER_IMAGE,
                     default_storage='500Mi', default_cpu=n_threads)
        window = '1e7' if pop == 'EUR' else '1e6'
        logger.info(f'Setting up {pop}...')
        chunk_size = int(5e6) if pop != 'EUR' else int(1e6)
        phenos_to_run = get_phenos_to_run(pop, limit=int(args.local_test),
                                          pilot=args.pilot,
                                          single_sex_only=args.single_sex_only, specific_phenos=args.phenos,
                                          skip_case_count_filter=args.skip_case_count_filter,
                                          first_round_phenos=args.run_first_round_phenos,
                                          sex_stratified=args.sex_stratified)
        logger.info(f'Got {len(phenos_to_run)} phenotypes...')
        if len(phenos_to_run) <= 20:
            logger.info(phenos_to_run)

        pheno_export_dir = f'{root}/pheno_export_data/{pop}'
        phenos_already_exported = {}
        if not args.overwrite_pheno_data and hl.hadoop_exists(pheno_export_dir):
            phenos_already_exported = {x['path'] for x in hl.hadoop_ls(pheno_export_dir)}
        pheno_exports = {}

        for pheno_key_dict in phenos_to_run:
            pheno_export_path = get_pheno_output_path(pheno_export_dir, pheno_key_dict, legacy=args.run_first_round_phenos)
            if not args.overwrite_pheno_data and pheno_export_path in phenos_already_exported:
                pheno_file = p.read_input(pheno_export_path)
            else:
                pheno_task = export_pheno(p, pheno_export_path, pheno_key_dict, 'ukbb_pan_ancestry',
                                          HAIL_DOCKER_IMAGE, additional_args=pop, n_threads=n_threads, proportion_single_sex=0)
                pheno_task.attributes.update({'pop': pop})
                pheno_file = pheno_task.out
            pheno_exports[stringify_pheno_key_dict(pheno_key_dict)] = pheno_file
        completed = Counter([isinstance(x, InputResourceFile) for x in pheno_exports.values()])
        logger.info(f'Exporting {completed[False]} phenos (already found {completed[True]})...')

        overwrite_null_models = args.create_null_models
        null_model_dir = f'{root}/null_glmm/{pop}'
        null_models_already_created = {}
        if not overwrite_null_models and hl.hadoop_exists(null_model_dir):
            null_models_already_created = {x['path'] for x in hl.hadoop_ls(null_model_dir)}
        null_models = {}

        for pheno_key_dict in phenos_to_run:
            null_glmm_root = get_pheno_output_path(null_model_dir, pheno_key_dict, '', legacy=args.run_first_round_phenos)
            model_file_path = f'{null_glmm_root}.rda'
            variance_ratio_file_path = f'{null_glmm_root}.{analysis_type}.varianceRatio.txt'

            if not overwrite_null_models and model_file_path in null_models_already_created and \
                    variance_ratio_file_path in null_models_already_created:
                model_file = p.read_input(model_file_path)
                variance_ratio_file = p.read_input(variance_ratio_file_path)
            else:
                if args.skip_any_null_models: continue
                fit_null_task = fit_null_glmm(p, null_glmm_root, pheno_exports[stringify_pheno_key_dict(pheno_key_dict)],
                                              pheno_key_dict['trait_type'], covariates,
                                              get_ukb_grm_plink_path(pop, iteration, window), SAIGE_DOCKER_IMAGE,
                                              inv_normalize=False, n_threads=n_threads, min_covariate_count=1,
                                              non_pre_emptible=args.non_pre_emptible, storage='100Gi')
                fit_null_task.attributes.update({'pop': pop})
                fit_null_task.attributes.update(copy.deepcopy(pheno_key_dict))
                model_file = fit_null_task.null_glmm.rda
                variance_ratio_file = fit_null_task.null_glmm[f'{analysis_type}.varianceRatio.txt']
            null_models[stringify_pheno_key_dict(pheno_key_dict)] = (model_file, variance_ratio_file)

        completed = Counter([isinstance(x[0], InputResourceFile) for x in null_models.values()])
        logger.info(f'Running {completed[False]} null models (already found {completed[True]})...')

        use_bgen = True
        vcf_dir = f'{root}/vcf/{pop}'
        test_extension = 'bgen' if use_bgen else 'vcf.gz'
        overwrite_vcfs = args.create_vcfs
        vcfs_already_created = {}
        if not overwrite_vcfs and hl.hadoop_exists(vcf_dir):
            vcfs_already_created = {x['path'] for x in hl.hadoop_ls(vcf_dir)}
        # logger.info(f'Found {len(vcfs_already_created)} VCFs in directory...')
        vcfs = {}
        for chromosome in chromosomes:
            chrom_length = chrom_lengths[chromosome]
            for start_pos in range(1, chrom_length, chunk_size):
                end_pos = chrom_length if start_pos + chunk_size > chrom_length else (start_pos + chunk_size)
                interval = f'{chromosome}:{start_pos}-{end_pos}'
                vcf_root = f'{vcf_dir}/variants_{chromosome}_{str(start_pos).zfill(9)}'
                if f'{vcf_root}.{test_extension}' in vcfs_already_created:
                    if use_bgen:
                        vcf_file = p.read_input_group(**{'bgen': f'{vcf_root}.bgen',
                                                         'bgen.bgi': f'{vcf_root}.bgen.bgi',
                                                         'sample': f'{vcf_root}.sample'})
                    else:
                        vcf_file = p.read_input_group(**{'vcf.gz': f'{vcf_root}.vcf.gz',
                                                         'vcf.gz.tbi': f'{vcf_root}.vcf.gz.tbi'})
                else:
                    vcf_task = extract_vcf_from_mt(p, vcf_root, HAIL_DOCKER_IMAGE, 'ukbb_pan_ancestry', adj=False,
                                                   additional_args=f'{chromosome},{pop}', input_dosage=True,
                                                   reference=reference, interval=interval, export_bgen=use_bgen,
                                                   n_threads=n_threads)
                    vcf_task.attributes['pop'] = pop
                    vcf_file = vcf_task.out
                vcfs[interval] = vcf_file
                if args.local_test:
                    break
            if args.local_test:
                break

        completed = Counter([type(x) == InputResourceFile for x in vcfs.values()])
        logger.info(f'Creating {completed[False]} VCFs (already found {completed[True]})...')

        result_dir = f'{root}/result/{pop}'
        overwrite_results = args.overwrite_results
        log_pvalue = True
        for i, pheno_key_dict in enumerate(phenos_to_run):
            if stringify_pheno_key_dict(pheno_key_dict) not in null_models: continue
            model_file, variance_ratio_file = null_models[stringify_pheno_key_dict(pheno_key_dict)]

            if not i % 10:
                n_jobs = dict(Counter(map(lambda x: x.name, p.select_jobs("")))).get("run_saige", 0)
                logger.info(f'Read {i} phenotypes ({n_jobs} new to run so far)...')

            pheno_results_dir = get_pheno_output_path(result_dir, pheno_key_dict, '', legacy=args.run_first_round_phenos)
            results_already_created = {}

            if not overwrite_results and not args.skip_saige and hl.hadoop_exists(pheno_results_dir):
                results_already_created = {x['path'] for x in hl.hadoop_ls(pheno_results_dir)}

            saige_tasks = []
            for chromosome in chromosomes:
                if args.skip_saige: break
                chrom_length = chrom_lengths[chromosome]
                for start_pos in range(1, chrom_length, chunk_size):
                    end_pos = chrom_length if start_pos + chunk_size > chrom_length else (start_pos + chunk_size)
                    interval = f'{chromosome}:{start_pos}-{end_pos}'
                    vcf_file = vcfs[interval]
                    results_path = get_results_prefix(pheno_results_dir, pheno_key_dict, chromosome, start_pos, legacy=args.run_first_round_phenos)
                    if overwrite_results or f'{results_path}.single_variant.txt' not in results_already_created:
                        samples_file = p.read_input(get_ukb_samples_file_path(pop, iteration))
                        saige_task = run_saige(p, results_path, model_file, variance_ratio_file, vcf_file, samples_file,
                                               SAIGE_DOCKER_IMAGE, trait_type=pheno_key_dict['trait_type'], use_bgen=use_bgen,
                                               chrom=chromosome, log_pvalue=log_pvalue)
                        saige_task.attributes.update({'interval': interval, 'pop': pop})
                        saige_task.attributes.update(copy.deepcopy(pheno_key_dict))
                        saige_tasks.append(saige_task)
                    if args.local_test:
                        break
                if args.local_test:
                    break

            res_tasks = []
            if overwrite_results or args.overwrite_hail_results or \
                    f'{pheno_results_dir}/variant_results.ht' not in results_already_created or \
                    not hl.hadoop_exists(f'{pheno_results_dir}/variant_results.ht/_SUCCESS'):
                null_glmm_root = get_pheno_output_path(null_model_dir, pheno_key_dict, f'.{analysis_type}.log',
                                                       legacy=args.run_first_round_phenos)

                prefix = get_results_prefix(pheno_results_dir, pheno_key_dict,
                                            f'{"chr" if reference == "GRCh38" else ""}{{chrom}}', 1,
                                            legacy=args.run_first_round_phenos)
                saige_log = f'{prefix}.{analysis_type}.log'

                load_task = load_results_into_hail(p, pheno_results_dir, pheno_key_dict,
                                                   saige_tasks, get_ukb_vep_path(), HAIL_DOCKER_IMAGE,
                                                   saige_log=saige_log, analysis_type=analysis_type,
                                                   n_threads=n_threads, null_glmm_log=null_glmm_root,
                                                   reference=reference, legacy_annotations=True,
                                                   log_pvalue=log_pvalue)
                load_task.attributes['pop'] = pop
                res_tasks.append(load_task)
                qq_export, qq_plot = qq_plot_results(p, pheno_results_dir, res_tasks, HAIL_DOCKER_IMAGE, QQ_DOCKER_IMAGE, n_threads=n_threads)
                qq_export.attributes.update({'pop': pop})
                qq_export.attributes.update(copy.deepcopy(pheno_key_dict))
                qq_plot.attributes.update({'pop': pop})
                qq_plot.attributes.update(copy.deepcopy(pheno_key_dict))

        logger.info(f'Setup took: {time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time))}')
        logger.info(f'Submitting: {get_tasks_from_pipeline(p)}')
        logger.info(f"Total size: {sum([len(x._pretty()) for x in p.select_jobs('')])}")
        p.run(dry_run=args.dry_run, wait=False, delete_scratch_on_exit=False)
        logger.info(f'Finished: {get_tasks_from_pipeline(p)}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite_pheno_data', help='Run single variant SAIGE', action='store_true')
    parser.add_argument('--single_sex_only', help='Run only single sex phenotypes (experimental)', action='store_true')
    parser.add_argument('--sex_stratified', help='Run these phenotypes in a sex-stratified fashion (experimental)', choices=(None, 'all', 'only'))
    parser.add_argument('--skip_any_null_models', help='Skip running SAIGE null models', action='store_true')
    parser.add_argument('--skip_saige', help='Skip running SAIGE tests', action='store_true')
    parser.add_argument('--create_null_models', help='Force creation of null models', action='store_true')
    parser.add_argument('--create_vcfs', help='Force creation of VCFs', action='store_true')
    parser.add_argument('--overwrite_results', help='Force run of SAIGE tests', action='store_true')
    parser.add_argument('--overwrite_hail_results', help='Force run of results loading', action='store_true')
    parser.add_argument('--local_test', help='Local test of pipeline', action='store_true')
    parser.add_argument('--non_pre_emptible', help='Local test of pipeline', action='store_true')
    parser.add_argument('--skip_case_count_filter', help='Skip running SAIGE tests', action='store_true')
    parser.add_argument('--phenos', help='Comma-separated list of trait_type-phenocode-pheno_sex-coding-modifier regexes '
                                         '(e.g. continuous-50-both_sexes--,icd10-E1.*,brain_mri-.* )')
    parser.add_argument('--pops', help='comma-searated list')
    parser.add_argument('--run_first_round_phenos', help='Run all phenotypes through pipeline (default: only pilot)', action='store_true')
    parser.add_argument('--pilot', help='Run all phenotypes through pipeline (default: only pilot)', action='store_true')
    parser.add_argument('--dry_run', help='Dry run only', action='store_true')
    args = parser.parse_args()

    main(args)
