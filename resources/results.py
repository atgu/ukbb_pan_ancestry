from .generic import *
from ukb_common.resources.results import *


def get_gene_intervals_path(reference: str = 'GRCh37'):
    return f'{public_bucket}/misc/gene_intervals_{reference}.ht'


def get_variant_results_path(pop: str, extension: str = 'mt'):
    return f'{bucket}/combined_results/results_{pop}.{extension}'


def get_variant_results_qc_path(extension: str = 'ht'):
    return f'{public_bucket}/sumstats_qc_analysis/full_variant_qc_metrics.{extension}'


def get_meta_analysis_results_path(extension: str = 'mt'):
    return f'{bucket}/combined_results/meta_analysis.{extension}'


def get_phenotype_results_qc_path(extension: str = 'ht'):
    return f'{bucket}/combined_results/full_phenotype_qc_metrics.{extension}'


def get_analysis_data_path(subdir: str, dataset: str, pop: str, extension: str = 'txt.bgz'):
    return f'{public_bucket}/sumstats_qc_analysis/{subdir}/{dataset}_{pop}.{extension}'


def get_results_timing_tsv_path(timing_type: str, pop: str = ''):
    check_timing_type(timing_type)

    pop = f'_{pop}' if timing_type == 'saige' else ''

    return f'{bucket}/results/misc/timings_{timing_type}{pop}.txt'


def get_results_timing_ht_path(timing_type: str):
    check_timing_type(timing_type)
    return f'{bucket}/results/misc/timings_{timing_type}.ht'


def get_heritability_txt_path(from_date: str = None):
    return f'{bucket}/results/misc/all_heritabilities{"_" + from_date if from_date else ""}.txt'