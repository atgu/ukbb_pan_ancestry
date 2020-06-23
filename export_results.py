#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 21:32:56 2020

Export flat file summary statistics from matrix tables

@author: nbaya
"""

import argparse
import hail as hl
from itertools import combinations
from time import time
from math import ceil
from ukbb_pan_ancestry.utils.results import load_final_sumstats_mt, get_meta_analysis_results_path, get_pheno_manifest_path

bucket = 'gs://ukb-diverse-pops'
public_bucket = 'gs://ukb-diverse-pops-public'
ldprune_dir = f'{bucket}/ld_prune'

def get_pheno_id(tb):
    pheno_id = (tb.trait_type+'-'+tb.phenocode+'-'+tb.pheno_sex+
                hl.if_else(hl.len(tb.coding)>0, '-'+tb.coding, '')+
                hl.if_else(hl.len(tb.modifier)>0, '-'+tb.modifier, '')
                ).replace(' ','_').replace('/','_')
    return pheno_id

def get_final_sumstats_mt_for_export():
    mt0 = load_final_sumstats_mt(filter_sumstats=False,
                                 filter_variants=False,
                                 separate_columns_by_pop=False,
                                 annotate_with_nearest_gene=False)
    mt0 = mt0.select_rows()
    return mt0

def export_results(num_pops, trait_types='all', batch_size=256, mt=None, 
                   export_path_str=None, skip_binary_eur=True):
    r'''
    `num_pops`: exact number of populations for which phenotype is defined
    `trait_types`: trait category (options: all, binary, quant)
    `batch_size`: batch size argument for export entries by col
    '''
    assert trait_types in {'all','quant','binary'}, "trait_types must be one of the following: {'all','quant','binary'}"
    print(f'\n\nExporting {trait_types} trait types for {num_pops} pops\n\n')
    if mt == None:
        mt0 = get_final_sumstats_mt_for_export()
    else:
        mt0 = mt
        
    meta_mt0 = hl.read_matrix_table(get_meta_analysis_results_path())
    
    mt0 = mt0.annotate_cols(pheno_id = get_pheno_id(tb=mt0))
    mt0 = mt0.annotate_rows(chr = mt0.locus.contig,
                            pos = mt0.locus.position,
                            ref = mt0.alleles[0],
                            alt = mt0.alleles[1])
    
    all_pops = ['AFR', 'AMR', 'CSA', 'EAS', 'EUR', 'MID']
    
    if trait_types == 'all':
        trait_types_to_run = ['continuous','biomarkers','categorical','phecode', 'icd10', 'prescriptions'] # list of which trait_type to run
    elif trait_types == 'quant':
        trait_types_to_run = ['continuous','biomarkers']
    elif trait_types == 'binary':
        trait_types_to_run = ['categorical','phecode', 'icd10', 'prescriptions']

    pop_sets = [set(i) for i in list(combinations(all_pops, num_pops))] # list of exact set of pops for which phenotype is defined
    
    
    # fields specific to each category of trait
    quant_meta_fields = ['AF_Allele2']
    quant_fields = ['AF_Allele2']
    
    binary_meta_fields = ['AF_Cases','AF_Controls']
    binary_fields = ['AF.Cases','AF.Controls']
    
    # dictionaries for renaming fields
    quant_meta_field_rename_dict = {'AF_Allele2':'af_meta',
                              'BETA':'beta_meta',
                              'SE':'se_meta',
                              'Pvalue':'pval_meta',
                              'Pvalue_het':'pval_heterogeneity'}
    quant_field_rename_dict = {'AF_Allele2':'af',
                         'BETA':'beta',
                         'SE':'se',
                         'Pvalue':'pval',
                         'low_confidence':'low_confidence'} # decided on this implementation to make later code cleaner
    
    binary_meta_field_rename_dict = {'BETA':'beta_meta',
                                     'SE':'se_meta',
                                     'Pvalue':'pval_meta',
                                     'AF_Cases':'af_cases_meta',
                                     'AF_Controls':'af_controls_meta',
                                     'Pvalue_het':'pval_heterogeneity'}
    binary_field_rename_dict = {'AF.Cases':'af_cases',
                                'AF.Controls':'af_controls',
                                'BETA':'beta',
                                'SE':'se',
                                'Pvalue':'pval',
                                'low_confidence':'low_confidence'} # decided on this implementation to make later code cleaner
    
    all_quant_trait_types = {'continuous','biomarkers'}
    all_binary_trait_types = {'categorical','phecode', 'icd10', 'prescriptions'}
    
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
            field_rename_dict = quant_field_rename_dict
        elif trait_category == 'binary':
            meta_fields = binary_meta_fields
            fields = binary_fields
            meta_field_rename_dict = binary_meta_field_rename_dict
            field_rename_dict = binary_field_rename_dict
    
        meta_fields += ['BETA','SE','Pvalue','Pvalue_het']
        fields += ['BETA','SE','Pvalue','low_confidence']
            
        for pop_set in pop_sets:    
            start = time()
            
            if (pop_set == {'EUR'} and trait_category == 'binary') and skip_binary_eur: # run EUR-only binary traits separately
                print('\nSkipping EUR-only binary traits\n')
                continue
            
            mt1 = mt0.filter_cols((hl.literal(trait_types).contains(mt0.trait_type))&
                                  (hl.set(mt0.pheno_data.pop)==hl.literal(pop_set)))
            
            col_ct = mt1.count_cols()
            if col_ct==0:
                print(f'\nSkipping {trait_types},{sorted(pop_set)}, no phenotypes found\n')
                continue
            
            pop_list = sorted(pop_set)
            
            annotate_dict = {}
            keyed_mt = meta_mt0[mt1.row_key,mt1.col_key]
            if len(pop_set)>1:
                for field in meta_fields: # NOTE: Meta-analysis columns go before per-population columns
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
            
            mt2 = mt1.annotate_entries(**annotate_dict)
            
            mt2 = mt2.filter_cols(mt2.coding != 'zekavat_20200409')
            mt2 = mt2.key_cols_by('pheno_id')
            mt2 = mt2.key_rows_by().drop('locus','alleles','summary_stats') # row fields that are no longer included: 'gene','annotation'
            
            batch_idx = 1        
            get_export_path = lambda batch_idx: f'{ldprune_dir}/export_results/{"" if export_path_str is None else f"{export_path_str}/"}{trait_category}/{"-".join(pop_list)}_batch{batch_idx}'
            print(mt2.describe())
            while hl.hadoop_is_dir(get_export_path(batch_idx)):
                batch_idx += 1
            print(f'\nExporting {col_ct} phenos to: {get_export_path(batch_idx)}\n')
            hl.experimental.export_entries_by_col(mt = mt2,
                                                  path = get_export_path(batch_idx),
                                                  bgzip = True,
                                                  batch_size = batch_size,
                                                  use_string_key_as_file_name = True,
                                                  header_json_in_file = False)
            end = time()
            print(f'\nExport complete for:\n{trait_types}\n{pop_list}\ntime: {round((end-start)/3600,2)} hrs')

def export_binary_eur(cluster_idx, num_clusters=10, batch_size = 256):
    r'''
    Export summary statistics for binary traits defined only for EUR. 
    Given the large number of such traits (4184), it makes sense to batch this 
    across `num_clusters` clusters for reduced wall time and robustness to mid-export errors.
    NOTE: `cluster_idx` is 1-indexed.
    '''
    mt0 = get_final_sumstats_mt_for_export()
    meta_mt0 = hl.read_matrix_table(get_meta_analysis_results_path())
    
    mt0 = mt0.annotate_cols(pheno_id = get_pheno_id(tb=mt0))
    mt0 = mt0.annotate_rows(chr = mt0.locus.contig,
                            pos = mt0.locus.position,
                            ref = mt0.alleles[0],
                            alt = mt0.alleles[1])
    
    trait_types_to_run = ['categorical','phecode', 'icd10', 'prescriptions'] # list of which trait_type to run
        
    # fields specific to each category of trait    
    meta_fields = ['AF_Cases','AF_Controls']
    fields = ['AF.Cases','AF.Controls']
    
    # dictionaries for renaming fields
    meta_field_rename_dict = {'BETA':'beta_meta',
                                     'SE':'se_meta',
                                     'Pvalue':'pval_meta',
                                     'AF_Cases':'af_cases_meta',
                                     'AF_Controls':'af_controls_meta',
                                     'Pvalue_het':'pval_heterogeneity'}
    field_rename_dict = {'AF.Cases':'af_cases',
                                'AF.Controls':'af_controls',
                                'BETA':'beta',
                                'SE':'se',
                                'Pvalue':'pval',
                                'low_confidence':'low_confidence'} # decided on this implementation to make later code cleaner
    
    all_binary_trait_types = {'categorical','phecode', 'icd10', 'prescriptions'}
    
    meta_fields += ['BETA','SE','Pvalue','Pvalue_het']
    fields += ['BETA','SE','Pvalue','low_confidence']
    
    trait_category = 'binary'        
    trait_types = all_binary_trait_types.intersection(trait_types_to_run) # get list of binary trait types to run
    pop_set = {'EUR'}
    start = time()
    
    mt1 = mt0.filter_cols((hl.literal(trait_types).contains(mt0.trait_type))&
                          (hl.set(mt0.pheno_data.pop)==hl.literal(pop_set)))
    
    pheno_id_list = mt1.pheno_id.collect()
    
    num_traits = len(pheno_id_list) # total number of traits to run
    
    traits_per_cluster = ceil(num_traits/num_clusters) # maximum traits to run per cluster
    
    cluster_pheno_id_list = pheno_id_list[(cluster_idx-1)*traits_per_cluster:cluster_idx*traits_per_cluster] # list of traits to run in current cluster
    
    print(len(cluster_pheno_id_list))
    
    mt1 = mt1.filter_cols(hl.literal(cluster_pheno_id_list).contains(mt1.pheno_id))
    
    pop_list = sorted(pop_set)
    
    annotate_dict = {}
    
    keyed_mt = meta_mt0[mt1.row_key,mt1.col_key]
    if len(pop_set)>1:
        for field in meta_fields: # NOTE: Meta-analysis columns go before per-population columns
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
    
    mt2 = mt1.annotate_entries(**annotate_dict)
    
    mt2 = mt2.filter_cols(mt2.coding != 'zekavat_20200409')
    mt2 = mt2.key_cols_by('pheno_id')
    mt2 = mt2.key_rows_by().drop('locus','alleles','summary_stats') # row fields that are no longer included: 'gene','annotation'
    print(mt2.describe())
    
    batch_idx = 1
    get_export_path = lambda batch_idx: f'{ldprune_dir}/release/{trait_category}/{"-".join(pop_list)}_batch{batch_idx}/subbatch{cluster_idx}'

    while hl.hadoop_is_dir(get_export_path(batch_idx)):
        batch_idx += 1
    print(f'\nExporting {len(cluster_pheno_id_list)} phenos to: {get_export_path(batch_idx)}\n')
    hl.experimental.export_entries_by_col(mt = mt2,
                                          path = get_export_path(batch_idx),
                                          bgzip = True,
                                          batch_size = batch_size,
                                          use_string_key_as_file_name = True,
                                          header_json_in_file = False)
    end = time()
    print(f'\nExport complete for:\n{trait_types}\n{pop_list}\ntime: {round((end-start)/3600,2)} hrs')
    
def export_subset(num_pops=None, phenocode=None):
    mt0 = get_final_sumstats_mt_for_export()
    if phenocode != None:
        print('\nFiltering to traits with phenocode: {phenocode}\n')
        mt0 = mt0.filter_cols(mt0.phenocode==phenocode)
    if num_pops is None:
        for num_pops in range(1,7):
                export_results(num_pops=num_pops, 
                               trait_types='all', 
                               batch_size=256, 
                               mt = mt0, 
                               export_path_str=phenocode)
    else:
        export_results(num_pops=num_pops, 
                       trait_types='all', 
                       batch_size=256, 
                       mt = mt0, 
                       export_path_str=phenocode)
                
    
def export_loo(batch_size=256, update=False):
    r'''
    For exporting p-values of meta-analysis of leave-one-out population sets
    '''
    meta_mt0 = hl.read_matrix_table(get_meta_analysis_results_path())
    
    meta_mt0 = meta_mt0.select_rows()

    meta_mt0 = meta_mt0.annotate_cols(pheno_id = get_pheno_id(tb=meta_mt0))
    
    meta_mt0 = meta_mt0.filter_cols(hl.len(meta_mt0.pheno_data.pop)==6)
    
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

    all_pops = sorted(['AFR', 'AMR', 'CSA', 'EAS', 'EUR', 'MID'])

    annotate_dict = {}
    '''
    pop_idx corresponds to the alphabetic ordering of the pops 
    entry with idx=0 is 6-pop meta-analysis, entry with idx=1 is 5-pop not-AFR 
    meta-analysis, idx=2 is 5-pop not-AMR, etc)
    '''
    for pop_idx, pop in enumerate(all_pops,1): 
        annotate_dict.update({f'pval_not_{pop}': meta_mt0.meta_analysis.Pvalue[pop_idx]})
    meta_mt1 = meta_mt0.annotate_entries(**annotate_dict)
    
    meta_mt1 = meta_mt1.key_cols_by('pheno_id')
    meta_mt1 = meta_mt1.key_rows_by().drop('locus','alleles','meta_analysis')
    
    meta_mt1.describe()
    
    batch_idx = 1
    get_export_path = lambda batch_idx: f'{ldprune_dir}/loo/sumstats/batch{batch_idx}'
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


def make_pheno_manifest(export=True):    
    mt0 = load_final_sumstats_mt(filter_sumstats=False,
                                 filter_variants=False,
                                 separate_columns_by_pop=False,
                                 annotate_with_nearest_gene=False)
    
    ht = mt0.cols()
    annotate_dict = {}
    
    annotate_dict.update({'pops': hl.delimit(ht.pheno_data.pop),
                          'num_pops': hl.len(ht.pheno_data.pop)})
     
    for field in ['n_cases','n_controls','heritability','lambda_gc']:
        for pop in ['AFR','AMR','CSA','EAS','EUR','MID']:
            new_field = field if field!='heritability' else 'saige_heritability' # new field name (only applicable to saige heritability)
            idx = ht.pheno_data.pop.index(pop)
            field_expr = ht.pheno_data[field]
            annotate_dict.update({f'{new_field}_{pop}': hl.if_else(hl.is_nan(idx),
                                                               hl.null(field_expr[0].dtype),
                                                               field_expr[idx])})
    annotate_dict.update({'filename': get_pheno_id(tb=ht)+'.tsv.bgz'})
    ht = ht.annotate(**annotate_dict)
    
    dropbox_manifest = hl.import_table(f'{ldprune_dir}/UKBB_Pan_Populations-Manifest_20200615-manifest_info.tsv',
                                       impute=True,
                                       key='File') # no need to filter table for duplicates because dropbox links are the same for updated phenos
    bgz = dropbox_manifest.filter(~dropbox_manifest.File.contains('.tbi'))
    bgz = bgz.rename({'File':'filename'})
    tbi = dropbox_manifest.filter(dropbox_manifest.File.contains('.tbi'))
    tbi = tbi.annotate(filename = tbi.File.replace('.tbi','')).key_by('filename')
        
    dropbox_annotate_dict = {}
    
    rename_dict = {'dbox link':'dropbox_link',
                   'size (bytes)':'size_in_bytes'}
    
    dropbox_annotate_dict.update({'filename_tabix':tbi[ht.filename].File})
    for field in ['dbox link','wget', 'size (bytes)','md5 hex']:
        for tb, suffix in [(bgz, ''), (tbi, '_tabix')]:
            dropbox_annotate_dict.update({(rename_dict[field] if field in rename_dict 
                                           else field.replace(' ','_')
                                           )+suffix:tb[ht.filename][field]})
    ht = ht.annotate(**dropbox_annotate_dict)
    ht = ht.drop('pheno_data')
    ht.describe()
    ht.show()
#    if export:
#        ht.export(get_pheno_manifest_path())
#    else:
#        return ht
    

if __name__=="__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--num-pops', type=int, default=None, help='number of defined pops (options: 6,5,...,2,1)')
    parser.add_argument('--trait-types', type=str, default='all', help='trait types to export (options: all, quant, binary)')
    parser.add_argument('--phenocode', type=str, default=None, help='phenocode to filter to for exporting results')
    parser.add_argument('--export-results', action='store_true')
    parser.add_argument('--export-binary_eur', action='store_true')
    parser.add_argument('--export-loo', action='store_true')
    parser.add_argument('--make-pheno-manifest', action='store_true')
    parser.add_argument('--export-updated-phenos', action='store_true')
    parser.add_argument('--cluster-idx',type=int, default=None, help='cluster index for splitting export of binary EUR traits')
    parser.add_argument('--batch_size', type=int, default=256, help='max number of phenotypes per batch for export_entries_by_col')
    args = parser.parse_args()

    if args.export_results:
        export_results(num_pops=args.num_pops,
                       trait_types=args.trait_types,
                       phenocode=args.phenocode)
    elif args.export_binary_eur:
        export_binary_eur(cluster_idx=args.cluster_idx)
    elif args.make_pheno_manifest:
        make_pheno_manifest()
    elif args.export_loo:
        export_loo(batch_size=args.batch_size)
    elif args.export_updated_phenos:
        export_updated_phenos(num_pops=args.num_pops)