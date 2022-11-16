#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 21:32:56 2020

Export flat file summary statistics from matrix tables

@author: nbaya, rahulg
"""

import argparse
from multiprocessing.sharedctypes import Value
import os, re
from copy import deepcopy
from sys import path
import hail as hl
import hailtop.batch as hb
from itertools import combinations
from time import time
from math import ceil

# hl.init(spark_conf={'spark.hadoop.fs.gs.requester.pays.mode': 'CUSTOM',
#                     'spark.hadoop.fs.gs.requester.pays.buckets': 'ukb-diverse-pops-public',
#                     'spark.hadoop.fs.gs.requester.pays.project.id': 'ukbb-diversepops-neale'})
# DISABLE WHEN CALLED
# hl.init(spark_conf={'spark.hadoop.fs.gs.requester.pays.mode': 'ENABLED',
#                      'spark.hadoop.fs.gs.requester.pays.project.id': 'ukbb-diversepops-neale'})

from ukbb_pan_ancestry.utils.results import load_final_sumstats_mt, load_meta_analysis_results
from ukbb_pan_ancestry.resources import POPS
from ukbb_pan_ancestry.resources.results import get_variant_results_path, get_pheno_manifest_path, get_h2_manifest_path, get_maximal_indepenedent_set_ht, get_h2_ht_path, get_h2_flat_file_path
from ukbb_pan_ancestry.resources.genotypes import get_filtered_mt
from ukbb_pan_ancestry.resources.phenotypes import get_ukb_pheno_mt_path
from ukbb_pan_ancestry.heritability.import_heritability import qc_to_flags
from ukbb_pan_ancestry.heritability.rhemc_pipeline import run_final_sink
from ukbb_pan_ancestry.get_timings_null_model import format_pheno_dir
from ukb_common.resources.generic import PHENO_KEY_FIELDS

bucket = 'gs://ukb-diverse-pops'
public_bucket = 'gs://ukb-diverse-pops-public'
ldprune_dir = f'{bucket}/ld_prune'
all_quant_trait_types = {'continuous','biomarkers'}
all_binary_trait_types = {'categorical','phecode', 'icd10', 'prescriptions'}
SEX = ['both_sexes','females','males']
path_max_indep_set = 'gs://ukb-diverse-pops/Cross-Pop-GWAS-comparison/Max_indep_set_phenos_h2QC_10_tiebreakCasenum_FINAL.ht'

# fields specific to each category of trait
quant_meta_fields = ['AF_Allele2']
quant_fields = ['AF_Allele2']

binary_meta_fields = ['AF_Cases','AF_Controls']
binary_fields = ['AF.Cases','AF.Controls']

# dictionaries for renaming fields
quant_meta_field_rename_dict = {'AF_Allele2': 'af_meta',
                                'BETA': 'beta_meta',
                                'SE': 'se_meta',
                                'Pvalue': 'pval_meta',
                                'Pvalue_het': 'pval_heterogeneity'}
quant_meta_hq_field_rename_dict = {'AF_Allele2': 'af_meta_hq',
                                   'BETA': 'beta_meta_hq',
                                   'SE': 'se_meta_hq',
                                   'Pvalue': 'pval_meta_hq',
                                   'Pvalue_het': 'pval_heterogeneity_hq'}
quant_field_rename_dict = {'AF_Allele2': 'af',
                           'BETA': 'beta',
                           'SE': 'se',
                           'Pvalue': 'pval',
                           'low_confidence': 'low_confidence'}

binary_meta_field_rename_dict = {'BETA': 'beta_meta',
                                 'SE': 'se_meta',
                                 'Pvalue': 'pval_meta',
                                 'AF_Cases': 'af_cases_meta',
                                 'AF_Controls': 'af_controls_meta',
                                 'Pvalue_het': 'pval_heterogeneity'}
binary_meta_hq_field_rename_dict = {'BETA': 'beta_meta_hq',
                                    'SE': 'se_meta_hq',
                                    'Pvalue': 'pval_meta_hq',
                                    'AF_Cases': 'af_cases_meta_hq',
                                    'AF_Controls': 'af_controls_meta_hq',
                                    'Pvalue_het': 'pval_heterogeneity_hq'}
binary_field_rename_dict = {'AF.Cases': 'af_cases',
                            'AF.Controls': 'af_controls',
                            'BETA': 'beta',
                            'SE': 'se',
                            'Pvalue': 'pval',
                            'low_confidence': 'low_confidence'}


def rename_dict_for_log10(dct, legacy_exp_p_values=False):
    dct_out = deepcopy(dct)
    if not legacy_exp_p_values:
        for k, v in dct_out.items():
            if re.search('^Pvalue.+', k) and re.search('^pval.+', v):
                dct_out.update({k: re.sub('^pval', 'neglog10_pval',v)})
    return dct_out


def get_pheno_id(tb):
    pheno_id = (tb.trait_type+'-'+tb.phenocode+'-'+tb.pheno_sex+
                hl.if_else(hl.len(tb.coding)>0, '-'+tb.coding, '')+
                hl.if_else(hl.len(tb.modifier)>0, '-'+tb.modifier, '')
                ).replace(' ','_').replace('/','_')
    return pheno_id


def get_final_sumstats_mt_for_export(exponentiate_p, custom_mt_path, legacy_exp_p_values):
    """ Updated to *not* filter by QC cutoffs.
    """
    mt0 = load_final_sumstats_mt(filter_sumstats=False,
                                 filter_variants=False,
                                 separate_columns_by_pop=False,
                                 annotate_with_nearest_gene=False,
                                 filter_pheno_h2_qc=False,
                                 exponentiate_p=exponentiate_p,
                                 legacy_exp_p_values=legacy_exp_p_values,
                                 custom_mt_path=custom_mt_path)
    mt0 = mt0.select_rows()
    return mt0


def export_results(num_pops, trait_types='all', batch_size=256, mt=None, 
                   export_path_str=None, skip_binary_eur=True, exponentiate_p=False,
                   legacy_exp_p_values=False,
                   suffix=None, skip_existing_folders=False, custom_mt_path=None):
    r'''
    `num_pops`: exact number of populations for which phenotype is defined
    `trait_types`: trait category (options: all, binary, quant)
    `batch_size`: batch size argument for export entries by col
    `suffix`: if not None, adds sumstats to a specified folder rather than 'export_results'
    '''
    assert trait_types in {'all','quant','binary'}, "trait_types must be one of the following: {'all','quant','binary'}"
    print(f'\n\nExporting {trait_types} trait types for {num_pops} pops\n\n')
    if mt == None:
        mt0 = get_final_sumstats_mt_for_export(exponentiate_p=exponentiate_p, custom_mt_path=custom_mt_path, legacy_exp_p_values=legacy_exp_p_values)
    else:
        mt0 = mt
        
    #meta_mt0 = hl.read_matrix_table(get_meta_analysis_results_path())
    meta_mt0 = load_meta_analysis_results(h2_filter='both', exponentiate_p=exponentiate_p, custom_path=os.path.dirname(custom_mt_path), legacy_exp_p_values=legacy_exp_p_values)
    
    mt0 = mt0.annotate_cols(pheno_id = get_pheno_id(tb=mt0))
    mt0 = mt0.annotate_rows(chr = mt0.locus.contig,
                            pos = mt0.locus.position,
                            ref = mt0.alleles[0],
                            alt = mt0.alleles[1])
    
    if trait_types == 'all':
        trait_types_to_run = ['continuous','biomarkers','categorical','phecode', 'icd10', 'prescriptions'] # list of which trait_type to run
    elif trait_types == 'quant':
        trait_types_to_run = ['continuous','biomarkers']
    elif trait_types == 'binary':
        trait_types_to_run = ['categorical','phecode', 'icd10', 'prescriptions']

    pop_sets = [set(i) for i in list(combinations(POPS, num_pops))] # list of exact set of pops for which phenotype is defined
    
    quant_trait_types = all_quant_trait_types.intersection(trait_types_to_run) # get list of quant trait types to run
    binary_trait_types = all_binary_trait_types.intersection(trait_types_to_run) # get list of binary trait types to run
    error_trait_types = set(trait_types_to_run).difference(quant_trait_types.union(binary_trait_types))
    assert len(error_trait_types)==0, f'ERROR: The following trait_types are invalid: {error_trait_types}'
        
    for trait_category, trait_types in [('binary', binary_trait_types), ('quant', quant_trait_types)]:
        if len(trait_types)==0: #if no traits in trait_types list
            continue

        print(f'{trait_category} trait types to run: {trait_types}')
        
        if trait_category == 'quant':
            meta_fields = quant_meta_fields
            fields = quant_fields
            meta_field_rename_dict = quant_meta_field_rename_dict
            meta_hq_field_rename_dict = quant_meta_hq_field_rename_dict
            field_rename_dict = quant_field_rename_dict
        elif trait_category == 'binary':
            meta_fields = binary_meta_fields
            fields = binary_fields
            meta_field_rename_dict = binary_meta_field_rename_dict
            meta_hq_field_rename_dict = binary_meta_hq_field_rename_dict
            field_rename_dict = binary_field_rename_dict

        meta_field_rename_dict = rename_dict_for_log10(meta_field_rename_dict, legacy_exp_p_values=legacy_exp_p_values)
        meta_hq_field_rename_dict = rename_dict_for_log10(meta_hq_field_rename_dict, legacy_exp_p_values=legacy_exp_p_values)
        field_rename_dict = rename_dict_for_log10(field_rename_dict, legacy_exp_p_values=legacy_exp_p_values)
    
        meta_fields += ['BETA','SE','Pvalue','Pvalue_het']
        fields += ['BETA','SE','Pvalue','low_confidence']
            
        for pop_set in pop_sets:
            
            get_export_path = lambda batch_idx: f'{ldprune_dir}/{"export_results" if suffix is None else suffix}/{"" if export_path_str is None else f"{export_path_str}/"}{trait_category}/{"-".join(pop_list)}_batch{batch_idx}'
            pop_list = sorted(pop_set)

            if skip_existing_folders:
                # We check if there are any folders with this set of ancestries; if so, skip
                path_to_export = os.path.dirname(get_export_path(1))
                if hl.hadoop_exists(path_to_export):
                    paths_found = [x['path'] for x in hl.hadoop_ls(path_to_export) if x['is_dir']]
                    anc_found = [re.sub('_batch[0-9]{1,}$','',os.path.basename(x)) for x in paths_found]
                    if "-".join(pop_list) in anc_found:
                        print(f'\nSkipping {"-".join(pop_list)} as its export folder was found\n')
                        continue
            
            start = time()
            
            if (pop_set == {'EUR'} and trait_category == 'binary') and skip_binary_eur: # run EUR-only binary traits separately
                print('\nSkipping EUR-only binary traits\n')
                continue
            
            mt1 = mt0.filter_cols((hl.literal(trait_types).contains(mt0.trait_type))&
                                  (hl.set(mt0.pheno_data.pop)==hl.literal(pop_set)))
            mt1 = mt1.annotate_cols(has_hq_meta_analysis = meta_mt0.cols()[mt1.col_key].has_hq_meta_analysis)
            
            col_ct = mt1.count_cols()
            if col_ct==0:
                print(f'\nSkipping {trait_types},{sorted(pop_set)}, no phenotypes found\n')
                continue

            # split into cols with meta analysis available and those without
            mt1_with_meta = mt1.filter_cols(hl.is_defined(mt1.has_hq_meta_analysis))
            mt1_no_meta = mt1.filter_cols(~hl.is_defined(mt1.has_hq_meta_analysis))
            n_no_meta = mt1_no_meta.count_cols()
            print(f'\nExporting {str(col_ct)} phenotypes...\n')
            print(f'\nNOTE: {str(n_no_meta)} phenotypes have no available meta analysis...\n')
            print(f'\nNOTE: {str(mt1_with_meta.count_cols())} phenotypes have an available meta analysis...\n')
            
            # we now split the mt into those with hq filtered meta analysis results and those without
            mt1_hq = mt1_with_meta.filter_cols(mt1_with_meta.has_hq_meta_analysis).drop('has_hq_meta_analysis')
            keyed_mt_hq_def = meta_mt0[mt1_hq.row_key,mt1_hq.col_key]

            mt1_hq_undef = mt1_with_meta.filter_cols(~mt1_with_meta.has_hq_meta_analysis).drop('has_hq_meta_analysis')
            keyed_mt_hq_undef = meta_mt0[mt1_hq_undef.row_key,mt1_hq_undef.col_key]


            def _shortcut_export_keyed(keyed_mt, mt1, use_hq, batch_idx):
                return _export_using_keyed_mt(keyed_mt, mt1=mt1, use_hq=use_hq, batch_idx=batch_idx,
                                              get_export_path=get_export_path,
                                              batch_size=batch_size, pop_set=pop_set,
                                              pop_list=pop_list, meta_fields=meta_fields, fields=fields,
                                              meta_field_rename_dict=meta_field_rename_dict,
                                              meta_hq_field_rename_dict=meta_hq_field_rename_dict,
                                              field_rename_dict=field_rename_dict)
            
            # export sumstats without hq columns
            if (n_no_meta > 0):
                batch_idx_nometa = _shortcut_export_keyed(None, mt1=mt1_no_meta, use_hq=None, batch_idx=1)
            else:
                batch_idx_nometa = 0

            # export sumstats with hq columns
            if (mt1_hq.count_cols() > 0):
                batch_idx_hq = _shortcut_export_keyed(keyed_mt_hq_def, mt1=mt1_hq, use_hq=True, batch_idx=batch_idx_nometa+1)
            else:
                batch_idx_hq = batch_idx_nometa
            
            # export sumstats without hq columns
            if (mt1_hq_undef.count_cols() > 0):
                _shortcut_export_keyed(keyed_mt_hq_undef, mt1=mt1_hq_undef, use_hq=False, batch_idx=batch_idx_hq+1)
            
            end = time()
            print(f'\nExport complete for:\n{trait_types}\n{pop_list}\ntime: {round((end-start)/3600,2)} hrs')


def export_binary_eur(cluster_idx, num_clusters=10, batch_size = 256, exponentiate_p=False,
                      suffix=None, custom_mt_path=None, legacy_exp_p_values=False):
    r'''
    Export summary statistics for binary traits defined only for EUR. 
    Given the large number of such traits (4184), it makes sense to batch this 
    across `num_clusters` clusters for reduced wall time and robustness to mid-export errors.
    NOTE: `cluster_idx` is 1-indexed.
    '''
    mt0 = get_final_sumstats_mt_for_export(exponentiate_p=exponentiate_p, custom_mt_path=custom_mt_path, legacy_exp_p_values=legacy_exp_p_values)
    #meta_mt0 = hl.read_matrix_table(get_meta_analysis_results_path())
    meta_mt0 = load_meta_analysis_results(h2_filter='both', exponentiate_p=exponentiate_p, custom_path=os.path.dirname(custom_mt_path), legacy_exp_p_values=legacy_exp_p_values)
    
    mt0 = mt0.annotate_cols(pheno_id = get_pheno_id(tb=mt0))
    mt0 = mt0.annotate_rows(chr = mt0.locus.contig,
                            pos = mt0.locus.position,
                            ref = mt0.alleles[0],
                            alt = mt0.alleles[1])
    
    trait_types_to_run = ['categorical','phecode', 'icd10', 'prescriptions'] # list of which trait_type to run
        
    # fields specific to each category of trait    
    meta_fields = binary_meta_fields + ['BETA','SE','Pvalue','Pvalue_het']
    fields = binary_fields + ['BETA','SE','Pvalue','low_confidence']

    trait_category = 'binary'        
    trait_types = all_binary_trait_types.intersection(trait_types_to_run) # get list of binary trait types to run
    pop_set = {'EUR'}
    start = time()
    
    mt1 = mt0.filter_cols((hl.literal(trait_types).contains(mt0.trait_type))&
                          (hl.set(mt0.pheno_data.pop)==hl.literal(pop_set)))
    mt1 = mt1.annotate_cols(has_hq_meta_analysis = meta_mt0.cols()[mt1.col_key].has_hq_meta_analysis)
    
    pheno_id_list = mt1.pheno_id.collect()
    
    num_traits = len(pheno_id_list) # total number of traits to run
    
    traits_per_cluster = ceil(num_traits/num_clusters) # maximum traits to run per cluster
    
    cluster_pheno_id_list = pheno_id_list[(cluster_idx-1)*traits_per_cluster:cluster_idx*traits_per_cluster] # list of traits to run in current cluster
    
    print(len(cluster_pheno_id_list))
    
    mt1 = mt1.filter_cols(hl.literal(cluster_pheno_id_list).contains(mt1.pheno_id))
    
    pop_list = sorted(pop_set)

    # we now split the mt into those with hq filtered meta analysis results and those without
    mt1_hq = mt1.filter_cols(mt1.has_hq_meta_analysis).drop('has_hq_meta_analysis')
    keyed_mt_hq_def = meta_mt0[mt1_hq.row_key,mt1_hq.col_key]

    mt1_hq_undef = mt1.filter_cols(~mt1.has_hq_meta_analysis).drop('has_hq_meta_analysis')
    keyed_mt_hq_undef = meta_mt0[mt1_hq_undef.row_key,mt1_hq_undef.col_key]

    get_export_path = lambda batch_idx: f'{ldprune_dir}/release{"" if suffix is None else "/"+suffix}/{trait_category}/{"-".join(pop_list)}_batch{batch_idx}/subbatch{cluster_idx}'
    
    
    meta_field_rename_dict = rename_dict_for_log10(binary_meta_field_rename_dict, legacy_exp_p_values=legacy_exp_p_values)
    meta_hq_field_rename_dict = rename_dict_for_log10(binary_meta_hq_field_rename_dict, legacy_exp_p_values=legacy_exp_p_values)
    field_rename_dict = rename_dict_for_log10(binary_field_rename_dict, legacy_exp_p_values=legacy_exp_p_values)

    def _shortcut_export_keyed(keyed_mt, mt1, use_hq, batch_idx):
        return _export_using_keyed_mt(keyed_mt, mt1=mt1, use_hq=use_hq, batch_idx=batch_idx,
                                      get_export_path=get_export_path,
                                      batch_size=batch_size, pop_set=pop_set,
                                      pop_list=pop_list, meta_fields=meta_fields, fields=fields,
                                      meta_field_rename_dict=meta_field_rename_dict,
                                      meta_hq_field_rename_dict=meta_hq_field_rename_dict,
                                      field_rename_dict=field_rename_dict)
    

    # export sumstats with hq columns
    if (mt1_hq.count_cols() > 0):
        batch_idx_hq = _shortcut_export_keyed(keyed_mt_hq_def, mt1=mt1_hq, use_hq=True, batch_idx=1)
    else:
        batch_idx_hq = 0
    
    # export sumstats without hq columns
    if (mt1_hq_undef.count_cols() > 0):
        _ = _shortcut_export_keyed(keyed_mt_hq_undef, mt1=mt1_hq_undef, use_hq=False, batch_idx=batch_idx_hq+1)

    end = time()
    print(f'\nExport complete for:\n{trait_types}\n{pop_list}\ntime: {round((end-start)/3600,2)} hrs')


def _export_using_keyed_mt(keyed_mt, mt1, use_hq, batch_idx, get_export_path,
                           batch_size, pop_set, pop_list, meta_fields, fields,
                           meta_field_rename_dict, meta_hq_field_rename_dict,
                           field_rename_dict):
    annotate_dict = {}
    if (keyed_mt is not None) and (len(pop_set)>1): # NOTE: Meta-analysis columns go before per-population columns
        if use_hq:
            for field in meta_fields:
                field_expr = keyed_mt.meta_analysis_hq[field][0]
                annotate_dict.update({f'{meta_hq_field_rename_dict[field]}': hl.if_else(hl.is_nan(field_expr),
                                                                                hl.str(field_expr),
                                                                                hl.format('%.3e', field_expr))})   
        for field in meta_fields:
            field_expr = keyed_mt.meta_analysis[field][0]
            annotate_dict.update({f'{meta_field_rename_dict[field]}': hl.if_else(hl.is_nan(field_expr),
                                                                    hl.str(field_expr),
                                                                    hl.format('%.3e', field_expr))})
    for field in fields:
        for pop_idx, pop in enumerate(pop_list):
            field_expr = mt1.summary_stats[field][pop_idx]
            annotate_dict.update({f'{field_rename_dict[field]}_{pop}': hl.if_else(hl.is_nan(field_expr),
                                                                    hl.str(field_expr),
                                                                    hl.str(field_expr) if field=='low_confidence' else hl.format('%.3e', field_expr))})
    
    new_ncol = mt1.count_cols()
    mt2 = mt1.annotate_entries(**annotate_dict)
    
    mt2 = mt2.filter_cols(mt2.coding != 'zekavat_20200409')
    mt2 = mt2.key_cols_by('pheno_id')
    # add 'chr','pos','ref','alt' to enforce ordering; should only matter if there are duplicate records
    mt2 = mt2.key_rows_by().drop('locus','alleles','summary_stats') # row fields that are no longer included: 'gene','annotation'
            
    print(mt2.describe())
    while hl.hadoop_is_dir(get_export_path(batch_idx)):
        batch_idx += 1
    print(f'\nExporting {new_ncol} phenos to: {get_export_path(batch_idx)}\n')
    hl.experimental.export_entries_by_col(mt = mt2,
                                          path = get_export_path(batch_idx),
                                          bgzip = True,
                                          batch_size = batch_size,
                                          use_string_key_as_file_name = True,
                                          header_json_in_file = False)
    return batch_idx


def load_phenotype_list(path):
    ht = hl.import_table(path).key_by(*PHENO_KEY_FIELDS)
    return ht


def export_subset(num_pops=None, phenocode=None, exponentiate_p=False, suffix=None,
                  skip_existing_folders=False, allow_binary_eur=False, legacy_exp_p_values=False,
                  export_specific_phenos=None, custom_mt_path=None):
    mt0 = get_final_sumstats_mt_for_export(exponentiate_p=exponentiate_p, custom_mt_path=custom_mt_path, legacy_exp_p_values=legacy_exp_p_values)
    if export_specific_phenos is not None:
        specific_ht = load_phenotype_list(export_specific_phenos)
        n_specific = specific_ht.count()
        mt0 = mt0.semi_join_cols(specific_ht)
        n_found = mt0.count_cols()
        print(f'Filtering to specific phenotypes in tsv. Of {str(n_specific)} phenotypes, {str(n_found)} were identified for exporting.')
    elif phenocode != None:
        print(f'\nFiltering to traits with phenocode: {phenocode}\n')
        mt0 = mt0.filter_cols(mt0.phenocode==phenocode)
    if num_pops is None:
        for num_pops in range(1,7):
                export_results(num_pops=num_pops, 
                               trait_types='all', 
                               batch_size=256, 
                               mt = mt0, 
                               export_path_str=phenocode,
                               exponentiate_p=exponentiate_p,
                               legacy_exp_p_values=legacy_exp_p_values,
                               suffix=suffix,
                               skip_existing_folders=skip_existing_folders,
                               skip_binary_eur = not allow_binary_eur,
                               custom_mt_path=custom_mt_path)
    else:
        export_results(num_pops=num_pops, 
                       trait_types='all', 
                       batch_size=256, 
                       mt = mt0, 
                       export_path_str=phenocode,
                       exponentiate_p=exponentiate_p,
                       legacy_exp_p_values=legacy_exp_p_values,
                       suffix=suffix,
                       skip_existing_folders=skip_existing_folders,
                       skip_binary_eur = not allow_binary_eur,
                       custom_mt_path=custom_mt_path)


def export_all_loo(batch_size=256, update=False, exponentiate_p=False, 
                   n_minimum_pops=3, suffix=None, h2_filter: bool=True,
                   export_specific_phenos=None, legacy_exp_p_values=False):
    """
    This function iterates through all phenotypes that have at least n_minimum_pops 
    and outputs loo meta-analysis results.

    NOTE not updated to use a custom mt
    """
    
    filter_string = 'pass' if h2_filter else 'none'
    meta_mt0 = load_meta_analysis_results(h2_filter=filter_string, exponentiate_p=exponentiate_p, legacy_exp_p_values=legacy_exp_p_values)   
    meta_mt0 = meta_mt0.select_rows()
    meta_mt0 = meta_mt0.annotate_cols(pheno_id = get_pheno_id(tb=meta_mt0))
    meta_mt0 = meta_mt0.filter_cols(hl.len(meta_mt0.pheno_data.pop)>=n_minimum_pops)

    if export_specific_phenos is not None:
        specific_ht = load_phenotype_list(export_specific_phenos)
        n_specific = specific_ht.count()
        meta_mt0 = meta_mt0.semi_join_cols(specific_ht)
        n_found = meta_mt0.count_cols()
        print(f'Filtering to specific phenotypes in tsv. Of {str(n_specific)} phenotypes, {str(n_found)} were identified for LOO exporting with {str(n_minimum_pops)} or more {"hq" if h2_filter else ""} pops.')
    
    if update:    
        current_dir = f'{ldprune_dir}/loo/sumstats/batch1' # directory of current results to update
        
        ss_list = hl.hadoop_ls(current_dir)
        pheno_id_list = [x['path'].replace('.tsv.bgz','').replace(f'{current_dir}/','') for x in ss_list if 'bgz' in x['path']]
    
        meta_mt0 = meta_mt0.filter_cols(~hl.literal(pheno_id_list).contains(meta_mt0.pheno_id))
                
    meta_mt0 = meta_mt0.annotate_rows(chr = meta_mt0.locus.contig,
                                      pos = meta_mt0.locus.position,
                                      SNP = (meta_mt0.locus.contig+':'+
                                             hl.str(meta_mt0.locus.position)+':'+
                                             meta_mt0.alleles[0]+':'+
                                             meta_mt0.alleles[1]))

    for num_pops in range(n_minimum_pops,7):
        pop_sets = [set(i) for i in list(combinations(sorted(POPS), num_pops))]
        for pop_set in pop_sets:
            pop_list = sorted(pop_set)
            meta_mt0_thisset = meta_mt0.filter_cols(meta_mt0.pheno_data.pop == hl.literal(pop_list))
            if meta_mt0_thisset.count_cols() > 0:
                export_loo(meta_mt0_thisset, batch_size=batch_size, pop_list=pop_list,
                           suffix=suffix, h2_filter=h2_filter)
            else:
                print(f'No {"hq " if h2_filter else ""}phenotypes found for {str(pop_list)}.')


def export_loo(meta_mt0, batch_size=256,
               pop_list=sorted(POPS), 
               suffix=None, h2_filter: bool=True):
    r'''
    For exporting p-values of meta-analysis of leave-one-out population sets.
    Now expects meta_mt0, which is generated in the export_all_loo function.
    To recreate the older version of this function, run export_all_loo(n_minimum_pops=6, h2_filter=False).
    '''
    meta_mt0 = meta_mt0.filter_cols(meta_mt0.pheno_data.pop == hl.literal(pop_list)) # for good measure
    annotate_dict = {}
    '''
    pop_idx corresponds to the alphabetic ordering of the pops. For 6 pops:
    entry with idx=0 is 6-pop meta-analysis, entry with idx=1 is 5-pop not-AFR 
    meta-analysis, idx=2 is 5-pop not-AMR, etc)
    '''
    for pop_idx, pop in enumerate(pop_list,1): 
        annotate_dict.update({f'pval_not_{pop}': meta_mt0.meta_analysis.Pvalue[pop_idx]})
    meta_mt1 = meta_mt0.annotate_entries(**annotate_dict)
    
    meta_mt1 = meta_mt1.key_cols_by('pheno_id')
    meta_mt1 = meta_mt1.key_rows_by().drop('locus','alleles','meta_analysis')
    
    batch_idx = 1
    get_export_path = lambda batch_idx: f'{ldprune_dir}/loo/sumstats/{"hq" if h2_filter else "all"}{"/"+suffix if suffix is not None else ""}/batch{batch_idx}'
    while hl.hadoop_is_dir(get_export_path(batch_idx)):
        batch_idx += 1
    print(f'\nExporting to: {get_export_path(batch_idx)}\n')
    print(meta_mt1.count_cols())
    hl.experimental.export_entries_by_col(mt = meta_mt1,
                                          path = get_export_path(batch_idx),
                                          bgzip = True,
                                          batch_size = batch_size,
                                          use_string_key_as_file_name = True,
                                          header_json_in_file = False)


def export_updated_phenos(num_pops=None):
    """
    NOTE not updated for new hq phenotypes
    """
    old_manifest = hl.import_table(get_pheno_manifest_path(),
                                   key=['trait_type','phenocode','pheno_sex',
                                        'coding','modifier'],
                                   impute=True)
    new_manifest = make_pheno_manifest(export=False)
    joined_manifest = old_manifest.join(new_manifest, how='outer')
    joined_manifest = joined_manifest.annotate(pheno_id = get_pheno_id(tb=joined_manifest))
    
    pheno_ids_new_phenos = joined_manifest.filter(hl.is_missing(joined_manifest.pops)).pheno_id.collect()
    pheno_ids_to_update = joined_manifest.filter(joined_manifest.pops!=joined_manifest.pops_1).pheno_id.collect()
    pheno_ids_new_phenos_str = '\n'.join(pheno_ids_new_phenos)
    pheno_ids_to_update_str = '\n'.join(pheno_ids_to_update)
    print(f'\n\nNew phenotypes to be exported:\n{pheno_ids_new_phenos_str}')
    print(f'\nUpdated phenotypes to be exported:\n{pheno_ids_to_update_str}')
    print(f'\n> Number of new phenotypes to be exported: {len(pheno_ids_new_phenos)}')
    print(f'\n> Number of phenotypes to be updated: {len(pheno_ids_to_update)}')
    print(f'\n> Total number of phenotypes to be exported: {len(pheno_ids_to_update+pheno_ids_new_phenos)}\n')
    
    # identify phenotypes that should be removed from previous results
    to_remove = joined_manifest.filter((hl.is_missing(joined_manifest.pops_1)))
    pheno_ids_to_remove = get_pheno_id(tb=to_remove).collect()
    pheno_ids_to_remove_str = '\n'.join(pheno_ids_to_remove)
    print(f'\nPhenotypes to remove:\n{pheno_ids_to_remove_str}')
    print(f'\n> Number of phenotypes to be removed: {len(pheno_ids_to_remove)}\n')
    
    # filtered to phenotypes that need to be updated (either completely new or a different set of populations)
    to_export = joined_manifest.filter((joined_manifest.pops!=joined_manifest.pops_1)|
                                       (hl.is_missing(joined_manifest.pops)))
    print(to_export.select('pops','pops_1').show(int(1e6)))
        
    mt0 = get_final_sumstats_mt_for_export()
    mt0 = mt0.filter_cols(hl.is_defined(to_export[mt0.col_key]))
    
    if num_pops == None:
        num_pops_set = set(to_export.num_pops_1.collect()) # get set of num_pops to run
        print(f'num_pops set: {num_pops_set}')
        for num_pops in num_pops_set:
            export_results(num_pops=num_pops, 
                           trait_types='all', 
                           batch_size=256, 
                           mt=mt0, 
                           export_path_str='update',
                           skip_binary_eur=False)
    else: # useful if parallelizing num_pops over multiple clusters
        export_results(num_pops=num_pops, 
                       trait_types='all', 
                       batch_size=256, 
                       mt=mt0, 
                       export_path_str='update',
                       skip_binary_eur=False)


def make_per_population_n_cases():
    """
    This function computes per-population and per-sex case count for each phenotype.
    """
    path_loc = f'{bucket}/combined_results/per_pop_per_sex_case_counts.ht'
    ht_saige = hl.read_matrix_table(get_variant_results_path('full')).cols()
    ht_saige = ht_saige.explode('pheno_data')
    ht_saige = ht_saige.key_by(*PHENO_KEY_FIELDS, pop = ht_saige.pheno_data.pop)

    if hl.hadoop_is_file(f'{path_loc}/_SUCCESS'):
        ht_per_pop = hl.read_table(path_loc)
    else:
        samples_keep = get_filtered_mt(chrom='22').cols()

        # number of samples
        # samples_keep.group_by(samples_keep.pop).aggregate(ct = hl.agg.count()).show()
        # +-------+--------+
        # | pop   |     ct |
        # +-------+--------+
        # | str   |  int64 |
        # +-------+--------+
        # | "AFR" |   6637 |
        # | "AMR" |    982 |
        # | "CSA" |   8876 |
        # | "EAS" |   2709 |
        # | "EUR" | 420542 |
        # | "MID" |   1599 |
        # +-------+--------+
        
        mt = hl.read_matrix_table(get_ukb_pheno_mt_path())
        mt = mt.key_rows_by(s = hl.str(mt.userId))
        mt = mt.semi_join_rows(samples_keep)

        # number of samples with phenotypes
        # ht_test = mt.rows()
        # ht_test.group_by(ht_test.pop).aggregate(ct = hl.agg.count()).show()
        # +-------+--------+
        # | pop   |     ct |
        # +-------+--------+
        # | str   |  int64 |
        # +-------+--------+
        # | "AFR" |   6636 |
        # | "AMR" |    980 |
        # | "CSA" |   8876 |
        # | "EAS" |   2709 |
        # | "EUR" | 420519 |
        # | "MID" |   1597 |
        # | NA    |     14 |
        # +-------+--------+
        # note this is almost identical to https://pan.ukbb.broadinstitute.org/docs/technical-overview with 12 fewer EUR and 2 fewer MID
        # except the 14 missing individuals are now "NA"; these seem like they were withdrawn as they are fully NA'd

        avoid_geq_0 = ['biomarkers', 'continuous', 'icd_first_occurrence']
        # if continuous, simply check if the entry is defined; if categorical, we also check if != 0 because we want to count cases
        mt = mt.annotate_cols(avoid_geq_0 = hl.literal(avoid_geq_0).contains(mt.trait_type))
        mt = mt.annotate_entries(**{f'count_{sex}': hl.is_defined(mt[sex]) & hl.if_else(mt.avoid_geq_0, True, mt[sex] != 0) for sex in SEX})
        ht_per_pop = mt.group_rows_by(mt.pop
                      ).aggregate(**{f'custom_n_cases_{sex}': hl.agg.count_where(mt[f'count_{sex}']) for sex in SEX}
                      ).entries(
                      ).drop('avoid_geq_0'
                      ).key_by(*PHENO_KEY_FIELDS, 'pop')
        ht_per_pop = ht_per_pop.persist()

        # now dealing with naming issues
        # - reverse "format_pheno_dir"
        # - add IRNT to biomarkers
        # - remove coding and modifier for specific continuous traits that are 3-keyed only in the saige cols table
        ht_per_pop_mapper = ht_per_pop.group_by(*PHENO_KEY_FIELDS).aggregate().key_by()
        saige_phenocode_mapper = ht_saige.group_by(*PHENO_KEY_FIELDS).aggregate().key_by()
        saige_phenocode_mapper = saige_phenocode_mapper.annotate(phenocode_new = format_pheno_dir(saige_phenocode_mapper.phenocode))
        saige_phenocode_mapper = saige_phenocode_mapper.key_by('phenocode_new')
        ht_per_pop_mapper = ht_per_pop_mapper.annotate(phenocode_updated=hl.case(missing_false=True)
                                                                           .when(hl.is_defined(saige_phenocode_mapper[ht_per_pop_mapper.phenocode]), 
                                                                                 saige_phenocode_mapper[ht_per_pop_mapper.phenocode].phenocode)
                                                                           .default(ht_per_pop_mapper.phenocode),
                                                modifier_updated=hl.case(missing_false=True)
                                                                   .when(ht_per_pop_mapper.trait_type == "biomarkers", "irnt")
                                                                   .when((ht_per_pop_mapper.trait_type == 'continuous') & (ht_per_pop_mapper.phenocode == 'random') & (ht_per_pop_mapper.modifier == 'random'), "")
                                                                   .default(ht_per_pop_mapper.modifier)).key_by(*PHENO_KEY_FIELDS)
        still_missing = saige_phenocode_mapper.key_by(*PHENO_KEY_FIELDS).anti_join(ht_per_pop_mapper.key_by('trait_type','phenocode_updated','pheno_sex','coding','modifier_updated'))
        # at this stage the only thing left should be continuous traits which have no coding or modifier in saige table
        if still_missing.filter((still_missing.trait_type != 'continuous') | (still_missing.modifier != "") | (still_missing.coding != "")).count() > 0:
            raise ValueError('ERROR: at this stage the only entries not mapping should be continuous with no modifier and no coding.')
        # now we just have to get rid of coding and modifier for these phenotypes
        ht_per_pop_mapper = ht_per_pop_mapper.key_by('trait_type', 'phenocode_updated', 'pheno_sex')
        still_missing = still_missing.key_by('trait_type', 'phenocode', 'pheno_sex')
        ht_per_pop_mapper_miss = ht_per_pop_mapper.semi_join(still_missing)
        # check that no keys are redundant
        if (ht_per_pop_mapper_miss.count() != ht_per_pop_mapper_miss.distinct().count()) | (still_missing.count() != still_missing.distinct().count()):
            raise ValueError('ERROR: in fixing names from phenotype table, some keys became duplicated.')
        # since we have 1:1 mapping from missing records to ht_per_pop_mapper, now just remove coding and modifier
        ht_per_pop_mapper = ht_per_pop_mapper.annotate(coding_updated=hl.case(missing_false=True)
                                                                        .when(hl.is_defined(still_missing[ht_per_pop_mapper.key]), "")
                                                                        .default(ht_per_pop_mapper.coding),
                                                       modifier_updated=hl.case(missing_false=True)
                                                                          .when(hl.is_defined(still_missing[ht_per_pop_mapper.key]), "")
                                                                          .default(ht_per_pop_mapper.modifier_updated))
        ht_per_pop_mapper_test = ht_per_pop_mapper.key_by('trait_type','phenocode_updated','pheno_sex','coding_updated','modifier_updated')
        # verify that there are no redundant keys
        if ht_per_pop_mapper_test.count() != ht_per_pop_mapper_test.distinct().count():
            raise ValueError('ERROR: in fixing names from phenotype table, some keys became duplicated.')

        originally_dupe_keys = ht_per_pop.group_by(**ht_per_pop.key).aggregate(ct= hl.agg.count())
        originally_dupe_keys = originally_dupe_keys.filter(originally_dupe_keys.ct > 1)
        
        ht_per_pop = ht_per_pop.key_by(*PHENO_KEY_FIELDS)
        ht_per_pop = ht_per_pop.annotate(fixed_keys = ht_per_pop_mapper.key_by(*PHENO_KEY_FIELDS)[ht_per_pop.key])
        ht_per_pop = ht_per_pop.key_by()
        ht_per_pop = ht_per_pop.annotate(phenocode = ht_per_pop.fixed_keys.phenocode_updated,
                                         coding = ht_per_pop.fixed_keys.coding_updated,
                                         modifier = ht_per_pop.fixed_keys.modifier_updated).key_by(*PHENO_KEY_FIELDS,'pop').drop('fixed_keys')
        
        newly_dupe_keys = ht_per_pop.group_by(**ht_per_pop.key).aggregate(ct= hl.agg.count())
        newly_dupe_keys = newly_dupe_keys.filter(newly_dupe_keys.ct > 1)
        new_dupes = originally_dupe_keys.anti_join(newly_dupe_keys.filter(hl.is_defined(newly_dupe_keys.pop)))

        if new_dupes.filter(hl.is_defined(new_dupes.pop)).count() > 0:
            raise ValueError('ERROR: there are new duplicate keys (including pop as a key, without undefined pop')
        if ht_saige.semi_join(newly_dupe_keys).count() > 0:
            raise ValueError('ERROR: duplicate keys should not be also found in Saige table as this could lead to problematic mapping.')
        if ht_saige.anti_join(ht_per_pop).count() > 0:
            raise ValueError('ERROR: not all records in saige table were found in the phenotype MatrixTable despite remapping efforts.')
        
        ht_per_pop = ht_per_pop.checkpoint(path_loc)


    ht_saige = ht_saige.select(n_cases_both_sexes = ht_saige.pheno_data.n_cases)
    ht = ht_saige.annotate(**{f'pre_qc_full_cohort_n_cases_{sex}': ht_per_pop[ht_saige.key][f'n_cases_{sex}'] for sex in SEX})
    ht = ht.annotate(**{f'recomputed_n_cases_{sex}': ht_per_pop[ht.key][f'custom_n_cases_{sex}'] for sex in SEX})
    
    # note that ht_saige had columns for n_cases_full_cohort_{sex}, but these are verifiably identical to those in ht_per_pop
    # run the following before dropping fields from ht_saige:
    # [ht.filter(ht[f'n_cases_full_cohort_{sex}'] != ht[f'pre_qc_full_cohort_n_cases_{sex}']).count() for sex in SEX]
    # >> [0, 0, 0]

    # below code can be used as a sanity check, since n_cases_both_sexes came from saige and recomputed_n_cases_both_sexes is from phenotype table
    # ht.filter(ht.n_cases_both_sexes != (ht.recomputed_n_cases_males + ht.recomputed_n_cases_females)).count()
    # >> 3121
    # ht.filter(ht.n_cases_both_sexes != (ht.recomputed_n_cases_both_sexes)).count()
    # >> 2116
    # ht.filter(~hl.is_defined(ht.pre_qc_full_cohort_n_cases_both_sexes)).count()
    # >> 0
    # ht.filter((ht.recomputed_n_cases_males == 0) & (ht.recomputed_n_cases_females == 0)).count()
    # >> 437
    # ht.filter((ht.recomputed_n_cases_males == 0) & (ht.recomputed_n_cases_females == 0) & ((ht.pre_qc_full_cohort_n_cases_males != 0) | (ht.pre_qc_full_cohort_n_cases_females != 0))).count()
    # >> 0
    
    return ht


def make_pheno_manifest(export=True, export_flattened_h2_table=False, web_version=False):    
    mt0 = load_final_sumstats_mt(filter_sumstats=False,
                                 filter_variants=False,
                                 separate_columns_by_pop=False,
                                 annotate_with_nearest_gene=False,
                                 filter_pheno_h2_qc=False)       
    
    ht = mt0.cols()
    ht = ht.annotate(phenotype_qc = hl.map(qc_to_flags, ht.pheno_data.heritability.qcflags))
    # DELETE BELOW TWO ROWS ONCE DONE VERIFYING THIS FUNCTION
    # ht = ht.rename({f'n_cases_full_cohort_{sex}': f'orig_n_cases_full_cohort_{sex}' for sex in SEX})
    # ht = ht.annotate(**{f'n_cases_full_cohort_{sex}': 1 for sex in SEX})
    #########
    ht = ht.annotate(**{f'n_cases_hq_cohort_{sex}': 1 for sex in SEX})
    annotate_dict = {}

    # pops passing QC are the same as pops in hq meta, so just keep pops passing QC
    #mt_meta = load_meta_analysis_results(h2_filter='pass')
    #ht_meta = mt_meta.cols()
    #ht_meta = ht_meta.annotate(pops_in_hq_meta = ht_meta.meta_analysis_data.pop[0])

    ht_max_indep = get_maximal_indepenedent_set_ht()

    annotate_dict.update({'pops': hl.delimit(ht.pheno_data.pop),
                          'num_pops': hl.len(ht.pheno_data.pop),
                          'pops_pass_qc': "",
                          'num_pops_pass_qc': 0})

    h2_fields = ['h2_observed','h2_observed_se','h2_liability','h2_liability_se','h2_z']
    for field in ['n_cases','n_controls',*h2_fields, 'lambda_gc', 'phenotype_qc']: # move saige h2 to the h2 table
        for pop in POPS:
            new_field = field if field!='saige_heritability' else 'saige_h2' # new field name (only applicable to saige heritability)
            prefix = ('sldsc_25bin_' if pop == 'EUR' else 'rhemc_25bin_50rv_') if new_field in h2_fields else ''
            new_field = prefix + new_field
            new_field = 'final_' + new_field if new_field in h2_fields else new_field
            idx = ht.pheno_data.pop.index(pop)
            if field == 'phenotype_qc':
                field_expr = ht[field]
            elif field in h2_fields:
                field_expr = ht.pheno_data.heritability.estimates.final[field]
            else:
                field_expr = ht.pheno_data[field]
            if (pop == 'AMR') and (field == 'phenotype_qc'):
                to_assn = hl.if_else(field_expr[idx] != 'GWAS_not_run', 'n_too_low', 'GWAS_not_run')
            else:
                to_assn = field_expr[idx]
            annotate_dict.update({f'{new_field}_{pop}': hl.if_else(hl.is_nan(idx),
                                                                   hl.null(field_expr[0].dtype),
                                                                   to_assn)})
    ht = ht.annotate(**annotate_dict)
    #ht = ht.annotate(pops_in_hq_meta = hl.delimit(ht_meta[ht.key].pops_in_hq_meta))
    #ht = ht.annotate(pops_in_hq_meta = hl.if_else(~hl.is_defined(ht.pops_in_hq_meta), "", ht.pops_in_hq_meta))
    ht = ht.annotate(pops_pass_qc_arr = hl.map(lambda y: y[0], hl.zip(ht.pheno_data.pop, ht.phenotype_qc).filter(lambda x: x[1] == 'PASS')))
    ht = ht.annotate(num_pops_pass_qc = hl.len(ht.pops_pass_qc_arr))
    ht = ht.annotate(pops_pass_qc = hl.delimit(ht.pops_pass_qc_arr))
    ht = ht.annotate(in_max_independent_set = ht_max_indep[ht.key].in_max_independent_set)
    ht = ht.annotate(in_max_independent_set = hl.if_else(hl.is_defined(ht.in_max_independent_set), ht.in_max_independent_set, False))
    ht = ht.annotate(filename = get_pheno_id(tb=ht)+'.tsv.bgz')
    ht = ht.annotate(filename_tabix = get_pheno_id(tb=ht)+'.tsv.bgz.tbi')
    if web_version:
        ht = ht.annotate(aws_link = 'https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/' + ht.filename)
        ht = ht.annotate(aws_link_tabix = 'https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files_tabix/' + ht.filename_tabix)
        ht = ht.annotate(wget = 'wget ' + ht.aws_link)
        ht = ht.annotate(wget_tabix = 'wget ' + ht.aws_link_tabix)
    else:
        ht = ht.annotate(aws_path = 's3://pan-ukb-us-east-1/sumstats_flat_files/' + ht.filename)
        ht = ht.annotate(aws_path_tabix = 's3://pan-ukb-us-east-1/sumstats_flat_files_tabix/' + ht.filename_tabix)

    # adding size/md5 for files and tabix files
    ht_size_md5 = hl.import_table(f'{bucket}/combined_results/2205_flat_file_info.tsv', impute=True, key='filename').rename({'md5':'md5_hex'})
    ht_size_md5_tbi = hl.import_table(f'{bucket}/combined_results/2205_tabix_file_info.tsv', impute=True, key='filename').rename({'md5':'md5_hex_tabix', 'size_in_bytes':'size_in_bytes_tabix'})
    ht = ht.annotate(**ht_size_md5[ht.filename])
    ht = ht.annotate(**ht_size_md5_tbi[ht.filename_tabix])


    # now create a table containing per-ancestry, per-sex case counts
    ht_counts = make_per_population_n_cases()
    ht_full_cohort = ht_counts.group_by(*PHENO_KEY_FIELDS
                             ).aggregate(n_cases_full_cohort_both_sexes = hl.agg.sum(ht_counts.n_cases_both_sexes), 
                                         n_cases_full_cohort_females = hl.agg.sum(ht_counts.recomputed_n_cases_females), 
                                         n_cases_full_cohort_males = hl.agg.sum(ht_counts.recomputed_n_cases_males))
    if ht.anti_join(ht_full_cohort).count() > 0:
        raise ValueError('ERROR: All keys in manifest should be found in the per-trait n_cases file.')
    ht = ht.annotate(**{f'n_cases_full_cohort_{sex}': ht_full_cohort[ht.key][f'n_cases_full_cohort_{sex}'] for sex in SEX})

    # do the same thing, now for the hq pops
    hq_pops_per_pheno = ht.select('pops_pass_qc_arr').explode('pops_pass_qc_arr').key_by(*PHENO_KEY_FIELDS, 'pops_pass_qc_arr')
    if hq_pops_per_pheno.anti_join(ht_counts).count() > 0:
        raise ValueError('ERROR: All hq population-trait should be found in the per-trait, per-population n_cases file.')
    ht_counts_f = ht_counts.semi_join(hq_pops_per_pheno)
    ht_counts_f.count()
    ht_hq_cohort = ht_counts_f.group_by(*PHENO_KEY_FIELDS
                             ).aggregate(n_cases_hq_cohort_both_sexes = hl.agg.sum(ht_counts_f.n_cases_both_sexes), 
                                         n_cases_hq_cohort_females = hl.agg.sum(ht_counts_f.recomputed_n_cases_females), 
                                         n_cases_hq_cohort_males = hl.agg.sum(ht_counts_f.recomputed_n_cases_males))
    ht = ht.annotate(**{f'n_cases_hq_cohort_{sex}': ht_hq_cohort[ht.key][f'n_cases_hq_cohort_{sex}'] for sex in SEX})

    # no longer using dropbox for distribution of sumstats
    # dropbox_manifest = hl.import_table(f'{ldprune_dir}/UKBB_Pan_Populations-Manifest_20200615-manifest_info.tsv',
    #                                    impute=True,
    #                                    key='File') # no need to filter table for duplicates because dropbox links are the same for updated phenos
    # bgz = dropbox_manifest.filter(~dropbox_manifest.File.contains('.tbi'))
    # bgz = bgz.rename({'File':'filename'})
    # tbi = dropbox_manifest.filter(dropbox_manifest.File.contains('.tbi'))
    # tbi = tbi.annotate(filename = tbi.File.replace('.tbi','')).key_by('filename')
        
    # dropbox_annotate_dict = {}
    
    # rename_dict = {'dbox link':'dropbox_link',
    #                'size (bytes)':'size_in_bytes'}
    
    # dropbox_annotate_dict.update({'filename_tabix':tbi[ht.filename].File})
    # for field in ['dbox link','wget', 'size (bytes)','md5 hex']:
    #     for tb, suffix in [(bgz, ''), (tbi, '_tabix')]:
    #         dropbox_annotate_dict.update({(rename_dict[field] if field in rename_dict 
    #                                        else field.replace(' ','_')
    #                                        )+suffix:tb[ht.filename][field]})
    # ht = ht.annotate(**dropbox_annotate_dict)
    ht = ht.drop('pheno_data', 'phenotype_qc', 'pops_pass_qc_arr')
    ht.describe()
    
    if export_flattened_h2_table:
        # now make h2 table
        ht_h2 = hl.import_table(get_h2_flat_file_path(), 
                                delimiter='\t', 
                                impute=True, 
                                key=PHENO_KEY_FIELDS)
        ht_h2 = ht_h2.rename({'ancestry':'pop'})
        ht_h2 = ht_h2.drop('phenotype_id')
        ht_h2 = ht_h2.key_by(*PHENO_KEY_FIELDS, 'pop')
        ht_annotate_saige = mt0.cols().explode('pheno_data')
        ht_annotate_saige = ht_annotate_saige.select(pop = ht_annotate_saige.pheno_data.pop,
                                                     saige_heritability = ht_annotate_saige.pheno_data.saige_heritability)
        ht_annotate_saige = ht_annotate_saige.key_by(*PHENO_KEY_FIELDS, 'pop')
        ht_h2 = ht_h2.annotate(**{'estimates.saige.h2': ht_annotate_saige[ht_h2.key].saige_heritability})
        ht_h2.describe()

    if export:
        #ht.export(get_pheno_manifest_path(web_version))
        ht.export(f'{bucket}/combined_results/220602_phenotype_manifest{"_web" if web_version else ""}.tsv.bgz')
        if export_flattened_h2_table:
            ht_h2.export(f'{bucket}/combined_results/220407_h2_manifest.tsv.bgz')
            #ht_h2.export(get_h2_manifest_path())
    else:
        return ht


def make_tabix(folder, sumstats_folder, suffix):
    if suffix is not None:
        folder = folder + suffix + '/'
        sumstats_folder = sumstats_folder + suffix + '/'
    
    backend = hb.ServiceBackend(billing_project='ukb_diverse_pops', bucket='ukb_diverse_pops')
    bserv = hb.Batch(name="tabix_pan_ancestry", backend=backend)
    if hl.hadoop_is_dir(sumstats_folder):
        items = [x['path'] for x in hl.hadoop_ls(sumstats_folder) if not x['is_dir']]
    else:
        raise ValueError('Sumstats folder does not exist.')
    if hl.hadoop_is_dir(folder):
        items_tabix = [x['path'] for x in hl.hadoop_ls(folder) if not x['is_dir']]
    else:
        items_tabix = []
    for sumstat in items:
        if folder + os.path.basename(sumstat) + '.tbi' not in items_tabix:
            j = bserv.new_job('tabix_' + os.path.basename(sumstat))
            j.image('gcr.io/ukbb-diversepops-neale/nbaya_tabix:latest')
            sumstat_file = bserv.read_input_group(**{'tsv.bgz': sumstat})
            command = f"""
                tabix -S1 -s1 -b2 -e2 {sumstat_file['tsv.bgz']}
                cp {sumstat_file['tsv.bgz']}.tbi {j.tbout}
            """
            j.command(command)
            bserv.write_output(j.tbout, folder + os.path.basename(sumstat) + '.tbi')
    bserv.run(verbose=False)


def get_md5_size_table(output_path, file_folder):
    """
    NOTE output_path should be the path with the filename WITHOUT extension.
    """
    backend = hb.ServiceBackend(billing_project='ukb_diverse_pops', remote_tmpdir='gs://ukbb-diverse-temp-7day/md5_size/')
    bserv = hb.Batch(name="md5_size_pan_ancestry", backend=backend)
    if hl.hadoop_is_dir(file_folder):
        items = [x['path'] for x in hl.hadoop_ls(file_folder) if not x['is_dir']]
    else:
        raise ValueError('File folder does not exist.')
    j_holder = []
    for item in items:
        j = bserv.new_job('size_md5_' + os.path.basename(item))
        j.image('gcr.io/ukbb-diversepops-neale/nbaya_tabix:latest')
        file = bserv.read_input(item)
        command = f"""
            thismd5=$(md5sum {file} | awk '{{print $1}}')
            thisbyte=$(wc -c {file} | awk '{{print $1}}')
            printf "filename\tmd5\tsize_in_bytes\n" > {j.tab_out}
            printf "{os.path.basename(item)}\t$thismd5\t$thisbyte\n" >> {j.tab_out}
        """
        j.command(command)
        j_holder.append(j)
    cat_script = bserv.read_input(f'gs://ukb-diverse-pops/rg-pcgc/code/concat_tables.py')
    _ = run_final_sink(bserv, j_holder, cat_script, 500, suffix='', output_file=os.path.basename(output_path), path_results=f'{os.path.dirname(output_path)}/')
    bserv.run(verbose=False)


if __name__=="__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--num-pops', type=int, default=None, help='number of defined pops (options: 6,5,...,2,1)')
    parser.add_argument('--trait-types', type=str, default='all', help='trait types to export (options: all, quant, binary)')
    parser.add_argument('--phenocode', type=str, default=None, help='phenocode to filter to for exporting results')

    parser.add_argument('--export-all-results', action='store_true', help='exports all results (except EUR binary) for all pops')
    parser.add_argument('--export-results', action='store_true')
    parser.add_argument('--export-binary_eur', action='store_true')
    parser.add_argument('--export-loo', action='store_true')
    parser.add_argument('--export-loo-minpops', type=int, default=3, help='smallest number of populations for which to output loo phenotypes')
    parser.add_argument('--export-loo-hq', action='store_true', help='if enabled along with --export-loo, will output only the hq loo results')
    parser.add_argument('--allow-binary-eur', action='store_true')
    parser.add_argument('--export-specific-phenos', type=str, default=None, help="path to a .tsv with specific 5-key phenotypes to export")

    parser.add_argument('--make-pheno-manifest', action='store_true')
    parser.add_argument('--export-h2-manifest', action='store_true')
    parser.add_argument('--export-web-manifest', action='store_true')

    parser.add_argument('--make-tabix', action='store_true')
    parser.add_argument('--tabix-folder', type=str, default='gs://ukb-diverse-pops/ld_prune/export_results/2112_final_results_tabix/')
    parser.add_argument('--sumstats-folder', type=str, default='gs://ukb-diverse-pops/ld_prune/export_results/2112_final_results/')
    
    parser.add_argument('--export-updated-phenos', action='store_true')
    parser.add_argument('--cluster-idx',type=int, default=None, help='cluster index for splitting export of binary EUR traits')
    parser.add_argument('--num-clusters',type=int, default=None, help='total number of clusters used in splitting export of binary EUR traits')
    parser.add_argument('--batch-size', type=int, default=256, help='max number of phenotypes per batch for export_entries_by_col')
    parser.add_argument('--exponentiate-p', action='store_true', help='enables regular scale p-values')
    parser.add_argument('--legacy-exp-p-values', action='store_true', help='If true, will revert to outputting exp(P). Default behavior is outputting -log10(P).')
    parser.add_argument('--suffix', type=str, default=None, help='if provided, will export to a folder specificed by suffix (added to default directory, so just give a folder name here')
    parser.add_argument('--custom-mt', type=str, default=None, help='if provided, will use this instead of the default pan-ukbb mt')
    parser.add_argument('--skip-existing-folders', action='store_true', help='for export_results and export_all_results, will skip a particular export if it exists (e.g., if quant/AFR_batch* exists, it is assumed the quant trait AFR export completed and it is skipped)')
    args = parser.parse_args()

    if args.export_results:
        export_results(num_pops=args.num_pops,
                       trait_types=args.trait_types,
                       batch_size=args.batch_size,
                       exponentiate_p=args.exponentiate_p,
                       legacy_exp_p_values=args.legacy_exp_p_values,
                       suffix=args.suffix, 
                       skip_existing_folders=args.skip_existing_folders,
                       custom_mt_path=args.custom_mt)
    elif args.export_all_results:
        # If phenocode is not provided, None will be provided below 
        # resulting in a full export across all pop combinations
        export_subset(exponentiate_p=args.exponentiate_p,
                      legacy_exp_p_values=args.legacy_exp_p_values,
                      phenocode=args.phenocode,
                      export_specific_phenos=args.export_specific_phenos,
                      suffix=args.suffix, num_pops=args.num_pops, allow_binary_eur=args.allow_binary_eur,
                      skip_existing_folders=args.skip_existing_folders,
                      custom_mt_path=args.custom_mt)
    elif args.export_binary_eur:
        export_binary_eur(batch_size=args.batch_size,
                          cluster_idx=args.cluster_idx,
                          num_clusters=args.num_clusters,
                          exponentiate_p=args.exponentiate_p,
                          legacy_exp_p_values=args.legacy_exp_p_values,
                          suffix=args.suffix,
                          custom_mt_path=args.custom_mt)
    elif args.make_pheno_manifest:
        make_pheno_manifest(export_flattened_h2_table=args.export_h2_manifest, web_version=args.export_web_manifest)
    elif args.export_loo:
        export_all_loo(batch_size=args.batch_size,
                       exponentiate_p=args.exponentiate_p,
                       legacy_exp_p_values=args.legacy_exp_p_values,
                       n_minimum_pops=args.export_loo_minpops,
                       suffix=args.suffix,
                       h2_filter=args.export_loo_hq,
                       export_specific_phenos=args.export_specific_phenos)
    elif args.export_updated_phenos:
        export_updated_phenos(num_pops=args.num_pops)
    elif args.make_tabix:
        make_tabix(args.tabix_folder, args.sumstats_folder, args.suffix)
