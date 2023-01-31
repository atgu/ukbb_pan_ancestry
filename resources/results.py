from .generic import *
from ukb_common.resources.results import *


def get_gene_intervals_path(reference: str = 'GRCh37'):
    return f'{public_bucket}/misc/gene_intervals_{reference}.ht'


def get_variant_results_path(pop: str, extension: str = 'mt'):
    if pop == 'full':
        return f'{public_bucket}/sumstats_release/results_{pop}.{extension}'
    else:
        return f'{bucket}/combined_results/results_{pop}.{extension}'


def get_variant_results_qc_path(extension: str = 'ht'):
    return f'{public_bucket_free}/sumstats_qc_analysis/full_variant_qc_metrics.{extension}'


def get_meta_analysis_results_path(filter_pheno_h2_qc: bool = True, extension: str = 'mt'):
    qc = 'h2_qc' if filter_pheno_h2_qc else 'raw'
    return f'{public_bucket}/sumstats_release/meta_analysis.{qc}.{extension}'


def get_phenotype_results_qc_path(extension: str = 'ht'):
    return f'{bucket}/combined_results/full_phenotype_qc_metrics.{extension}'


def get_h2_flat_file_path():
    return f'{public_bucket_free}/h2/h2_estimates_all_flat_211101.tsv'


def get_h2_ht_path():
    return f'{public_bucket_free}/h2/h2_estimates_all.ht'


def get_maximal_independent_set_path():
    return f'{public_bucket_free}/qc/Max_indep_set_phenos_h2QC_10_tiebreakCasenum.ht'


def get_maximal_indepenedent_set_ht():
    return hl.read_table(get_maximal_independent_set_path()).annotate(in_max_independent_set=True)


def get_analysis_data_path(subdir: str, dataset: str, pop: str, extension: str = 'txt.bgz'):
    return f'{public_bucket_free}/sumstats_qc_analysis/{subdir}/{dataset}_{pop}.{extension}'


def get_final_lambdas_path():
    return get_analysis_data_path('lambda', 'lambdas', 'full', 'ht')


def get_results_timing_tsv_path(timing_type: str, pop: str = ''):
    check_timing_type(timing_type)

    pop = f'_{pop}' if timing_type == 'saige' else ''

    return f'{bucket}/results/misc/timings_{timing_type}{pop}.txt'


def get_results_timing_ht_path(timing_type: str):
    check_timing_type(timing_type)
    return f'{bucket}/results/misc/timings_{timing_type}.ht'


def get_heritability_txt_path(from_date: str = None):
    return f'{bucket}/results/misc/all_heritabilities{"_" + from_date if from_date else ""}.txt'


def get_pheno_manifest_path(web_version=False):
    return f'{public_bucket}/sumstats_release/phenotype_manifest{"_web" if web_version else ""}.tsv.bgz'


def get_h2_manifest_path():
    return f'{public_bucket}/sumstats_release/h2_manifest.tsv.bgz'


def get_clumping_results_path(pop: str = 'full', high_quality: bool = False, 
                              not_pop: bool = True, max_pops: bool = False):
    
    """
    Clumping results available for only high quality variants and high quality phenotypes (single pops, leave-one-out meta-analyses, and all-pop meta-analyses).
    
    :param pop: Input pop for single pop or leave-one-out results or "full" for all results in one MT
    :param high_quality: High quality variants only
    :param not_pop: Leave-one-out meta-analyses
    :param max_pops: All-pop meta-analyses
    :param hq_phenos: Leave-one-out meta-analyses (summstats available only for hq phenos)
    
    """
        
    mt_name = 'max_pops' if max_pops else (f'{"not_" if not_pop else ""}{pop}' if pop in POPS else 'full_clump_results')
    
    if high_quality:
        return f'{ldprune_dir}/clump_results_high_quality_22115/{mt_name}.mt' 
    else:
        return f'{ldprune_dir}/clump_results/{mt_name}.mt'

def get_prs_mt_path(high_quality: bool = True):
    return f'{bucket}/prs/all_combos_prs{"" if high_quality else "_raw"}.mt'


def get_clump_sumstats_bm_path(high_quality: bool = True):
    return f'{public_bucket}/misc/prs/clumped_sumstats{"" if high_quality else "_raw"}.bm'


def get_clump_sumstats_col_ht_path(high_quality: bool = True):
    return f'{public_bucket}/misc/prs/clumped_sumstats{"" if high_quality else "_raw"}.cols.ht'


genotype_bm_path = f'{temp_bucket}/prs/genotypes.bm'
genotype_samples_ht_path = f'{temp_bucket}/prs/genotype_samples.ht'


def get_prs_bm_path(high_quality: bool = True):
    return f'{temp_bucket}/prs/prs{"" if high_quality else "_raw"}.bm'


def get_prs_assess_ht_path(high_quality: bool = True):
    return f'{bucket}/prs/assess_prs{"" if high_quality else "_raw"}.ht'
