from .generic import *
from ukb_common.resources.results import *

ldprune_dir = f'{bucket}/ld_prune'

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
    return f'{public_bucket_free}/h2/h2_estimates_all_flat_221123.tsv'


def get_h2_ht_path():
    return f'{public_bucket_free}/h2/h2_estimates_all.ht'


def get_maximal_independent_set_path():
    return f'{public_bucket_free}/qc/Max_indep_set_phenos_h2QC_10_tiebreakCasenum_FINAL.ht'


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
    
    :param pop: Input pop for single pop or leave-one-out results or "full" for all results in one MT (single pop summstats available for all phenos; clumping results available for only hq phenos)
    :param high_quality: High quality variants only
    :param not_pop: Leave-one-out meta-analyses (summstats available for only hq phenos; clumping results available for only hq phenos)
    :param max_pops: All-pop meta-analyses (summstats available for all phenos; clumping results available for only hq phenos)
    
    """
        
    mt_name = ("max_pops" if max_pops else ("not_" if not_pop else ""))+(f'{pop}' if pop in POPS else ('full' if pop == "full" and not max_pops else ""))
    
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


def get_distance_clumping_results_path(pop: str, radius: int, merged: bool = True, extension: str = "ht"):
    merged = "merged" if merged else "not_merged"
    return f"{bucket}/misc/distance_clump_results.{pop}.r{radius}.{merged}.{extension}"


genotype_bm_path = f'{temp_bucket}/prs/genotypes.bm'
genotype_samples_ht_path = f'{temp_bucket}/prs/genotype_samples.ht'


def get_prs_bm_path(high_quality: bool = True):
    return f'{temp_bucket}/prs/prs{"" if high_quality else "_raw"}.bm'


def get_prs_assess_ht_path(high_quality: bool = True):
    return f'{bucket}/prs/assess_prs{"" if high_quality else "_raw"}.ht'

otg_release = "22.09"


def get_ukb_pheno_efo_mapping_path(release: str = otg_release, extension: str = "ht"):
    return f'{bucket}/misc/otg/{release}/ukb_pheno_efo_mapping.{extension}'


def get_munged_otg_v2d_path(release: str = otg_release, extension: str = "ht"):
    return f'{bucket}/misc/otg/{release}/munged_otg_v2d.{extension}'


def get_known_ukbb_loci_path(pop: str, release: str = otg_release, extension: str = "ht"):
    return f'{bucket}/misc/otg/{release}/known_ukbb_loci.{pop}.{extension}'
