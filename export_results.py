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
from diverse_pops.utils.results import get_analysis_data_path, load_final_sumstats_mt, get_meta_analysis_results_path


bucket = 'gs://ukb-diverse-pops'
public_bucket = 'gs://ukb-diverse-pops-public'
ldprune_dir = f'{bucket}/ld_prune'


def export_results(num_pops, trait_types='all', batch_size=256, phenocode=None):
    r'''
    `num_pops`: exact number of populations for which phenotype is defined
    `trait_types`: trait category (options: all, binary, quant)
    `batch_size`: batch size argument for export entries by col
    '''
    assert trait_types in {'all','quant','binary'}, "trait_types must be one of the following: {'all','quant','binary'}"
    hl.init(default_reference='GRCh38', log='/tmp/export_entries_by_col.log')
    print(f'\n\nExporting {trait_types} trait types for {num_pops} pops\n\n')
    mt0 = load_final_sumstats_mt(filter_sumstats=False,
                                 filter_variants=False,
                                 separate_columns_by_pop=False,
                                 annotate_with_nearest_gene=False)
    mt0 = mt0.select_rows()
    
    if phenocode != None:
        print('\nFiltering to traits with phenocode: {phenocode}\n')
        mt0 = mt0.filter_cols(mt0.phenocode==phenocode)
        
    meta_mt0 = hl.read_matrix_table(get_meta_analysis_results_path())
    
    mt0 = mt0.annotate_cols(pheno_id = (mt0.trait_type+'-'+
                            mt0.phenocode+'-'+
                            mt0.pheno_sex+
                            hl.if_else(hl.len(mt0.coding)>0, '-'+mt0.coding, '')+
                            hl.if_else(hl.len(mt0.modifier)>0, '-'+mt0.modifier, '')
                            ).replace(' ','_').replace('/','_'))
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
            
            if pop_set == {'EUR'} and trait_category == 'binary': # run EUR-only binary traits separately
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
            # TODO: Filter variants with NA in field
            keyed_mt = meta_mt0[mt1.row_key,mt1.col_key]
            if len(pop_set)>1:
                for field in meta_fields: # NOTE: Meta-analysis columns go before per-population columns
            #        annotate_dict.update({f'{meta_field_rename_dict[field]}': hl.float64(hl.format('%.3e', keyed_mt.meta_analysis[field][0]))})
                    field_expr = keyed_mt.meta_analysis[field][0]
                    annotate_dict.update({f'{meta_field_rename_dict[field]}': hl.if_else(hl.is_nan(field_expr),
                                                                              hl.str(field_expr),
                                                                              hl.format('%.3e', field_expr))})
        
            for field in fields:
                for pop_idx, pop in enumerate(pop_list):
        #            annotate_dict.update({f'{field_rename_dict[field]}_{pops[pop_idx]}': hl.format('%.3e', mt1.summary_stats[field][pop_idx])})
                    field_expr = mt1.summary_stats[field][pop_idx]
                    annotate_dict.update({f'{field_rename_dict[field]}_{pop}': hl.if_else(hl.is_nan(field_expr),
                                                                               hl.str(field_expr),
                                                                               hl.str(field_expr) if field=='low_confidence' else hl.format('%.3e', field_expr))})
            
            mt2 = mt1.annotate_entries(**annotate_dict)
            
            mt2 = mt2.filter_cols(mt2.coding != 'zekavat_20200409')
            mt2 = mt2.key_cols_by('pheno_id')
            mt2 = mt2.key_rows_by().drop('locus','alleles','summary_stats') # row fields that are no longer included: 'gene','annotation'
            
            batch_idx = 1        
            get_export_path = lambda batch_idx: f'{ldprune_dir}/release/{trait_category}/{"" if phenocode is None else f"{phenocode}/"}{"-".join(pop_list)}_batch{batch_idx}'
#            if hl.hadoop_is_dir(get_export_path(batch_idx)):
#                print(f'Skipping because path exists: {get_export_path(batch_idx)}')
#                continue        
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
    
    hl.init(default_reference='GRCh38', log='/tmp/export_entries_by_col.log')
    mt0 = load_final_sumstats_mt(filter_sumstats=False,
                                 filter_variants=False,
                                 separate_columns_by_pop=False,
                                 annotate_with_nearest_gene=False)
    mt0 = mt0.select_rows()
    meta_mt0 = hl.read_matrix_table(get_meta_analysis_results_path())
    
    
    mt0 = mt0.annotate_cols(pheno_id = (mt0.trait_type+'-'+
                            mt0.phenocode+'-'+
                            mt0.pheno_sex+
                            hl.if_else(hl.len(mt0.coding)>0, '-'+mt0.coding, '')+
                            hl.if_else(hl.len(mt0.modifier)>0, '-'+mt0.modifier, '')
                            ).replace(' ','_').replace('/','_'))
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
    #        annotate_dict.update({f'{meta_field_rename_dict[field]}': hl.float64(hl.format('%.3e', keyed_mt.meta_analysis[field][0]))})
            field_expr = keyed_mt.meta_analysis[field][0]
            annotate_dict.update({f'{meta_field_rename_dict[field]}': hl.if_else(hl.is_nan(field_expr),
                                                                      hl.str(field_expr),
                                                                      hl.format('%.3e', field_expr))})

    for field in fields:
        for pop_idx, pop in enumerate(pop_list):
#            annotate_dict.update({f'{field_rename_dict[field]}_{pops[pop_idx]}': hl.format('%.3e', mt1.summary_stats[field][pop_idx])})
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
    
def export_loo(batch_size=256):
    r'''
    For exporting p-values of meta-analysis of leave-one-out population sets
    '''
    meta_mt0 = hl.read_matrix_table(get_meta_analysis_results_path())
    
    meta_mt0 = meta_mt0.filter_cols(hl.len(meta_mt0.pheno_data.pop)==6)
    
    meta_mt0 = meta_mt0.annotate_cols(pheno_id = (meta_mt0.trait_type+'-'+
                                                  meta_mt0.phenocode+'-'+
                                                  meta_mt0.pheno_sex+
                                                  hl.if_else(hl.len(meta_mt0.coding)>0, '-'+meta_mt0.coding, '')+
                                                  hl.if_else(hl.len(meta_mt0.modifier)>0, '-'+meta_mt0.modifier, '')
                                                  ).replace(' ','_').replace('/','_'))
                
    meta_mt0 = meta_mt0.annotate_rows(SNP = (meta_mt0.locus.contig+':'+
                                             hl.str(meta_mt0.locus.position)+':'+
                                             meta_mt0.alleles[0]+':'+
                                             meta_mt0.alleles[1]))

    all_pops = sorted(['AFR', 'AMR', 'CSA', 'EAS', 'EUR', 'MID'])
    
    annotate_dict = {}
    for pop_idx, pop in enumerate(all_pops,1): # pop idx corresponds to the alphabetic ordering of the pops (entry with idx=0 is 6-pop meta-analysis)
        annotate_dict.update({f'pval_not_{pop}': meta_mt0.meta_analysis.Pvalue[pop_idx]})
    meta_mt1 = meta_mt0.annotate_entries(**annotate_dict)
    
    meta_mt1 = meta_mt1.key_cols_by('pheno_id')
    meta_mt1 = meta_mt1.key_rows_by().drop('locus','alleles','gene','annotation','meta_analysis')
    
    print(meta_mt1.describe())
    
    batch_idx = 1
    get_export_path = lambda batch_idx: f'{ldprune_dir}/loo/sumstats/batch{batch_idx}'
    while hl.hadoop_is_dir(get_export_path(batch_idx)):
        batch_idx += 1
    print(f'\nExporting to: {get_export_path(batch_idx)}\n')
    hl.experimental.export_entries_by_col(mt = meta_mt1,
                                          path = get_export_path(batch_idx),
                                          bgzip = True,
                                          batch_size = batch_size,
                                          use_string_key_as_file_name = True,
                                          header_json_in_file = False)
def make_pheno_manifest():    
    mt0 = load_final_sumstats_mt(filter_sumstats=False,
                                 filter_variants=False,
                                 separate_columns_by_pop=False,
                                 annotate_with_nearest_gene=False)
    
#    ht = hl.read_table(get_analysis_data_path('lambda', 'lambdas', 'full', 'ht')) # using this table results in 7256 phenos instead of 720
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
    annotate_dict.update({'filename': (ht.trait_type+'-'+
                                          ht.phenocode+'-'+
                                          ht.pheno_sex+
                                          hl.if_else(hl.len(ht.coding)>0, '-'+ht.coding, '')+
                                          hl.if_else(hl.len(ht.modifier)>0, '-'+ht.modifier, '')
                                          ).replace(' ','_').replace('/','_')+
                                          '.tsv.bgz'})
    ht = ht.annotate(**annotate_dict)
#    dropbox_dir= f'{bucket}/sumstats_flat_files'
#    ht = ht.annotate('file_location' = dropbox_dir+'/'+})
    ht = ht.drop('pheno_data','pheno_indices')
    ht.export(f'{bucket}/combined_results/phenotype_manifest.tsv.bgz')
    

if __name__=="__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--num-pops', type=int, default=0, help='number of defined pops (options: 6,5,...,2,1)')
    parser.add_argument('--trait-types', type=str, default='all', help='trait types to export (options: all, quant, binary)')
    parser.add_argument('--phenocode', type=str, default=None, help='phenocode to filter to for exporting results')
    parser.add_argument('--export-results', action='store_true')
    parser.add_argument('--export-binary_eur', action='store_true')
    parser.add_argument('--export-loo', action='store_true')
    parser.add_argument('--make-pheno-manifest', action='store_true')
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