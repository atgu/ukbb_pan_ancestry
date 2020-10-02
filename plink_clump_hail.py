#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 15:09:04 2020

Hail script for clumping GWAS results with PLINK

@author: nbaya
"""

import argparse
import hail as hl
import sys
import ukb_common
from ukbb_pan_ancestry import get_clumping_results_path #get_pheno_manifest_path
#from ukb_common import mwzj_hts_by_tree

bucket = 'gs://ukb-diverse-pops'
ldprune_dir = f'{bucket}/ld_prune'

all_pops = ['AFR', 'AMR', 'CSA', 'EAS', 'EUR', 'MID']

def ht_to_tsv(args):
    r'''
    Convert Hail table of variant results to a tsv 
    '''
    add_args = {}
    if args.n_threads is not None:
        add_args['master'] = f'local[{args.n_threads}]'
    hl.init(default_reference='GRCh38', log='/tmp/ht_to_tsv.log', **add_args)
    ht = hl.read_table(args.input_file)
    print(ht.describe())
    ht.export(args.output_file)

    
def tsv_to_ht(args):
    r'''
    Convert tsv to a Hail table
    '''
    print(sys.version)
    add_args = {}
    if args.n_threads is not None:
        add_args['master'] = f'local[{args.n_threads}]'
    hl.init(default_reference='GRCh38', log='/tmp/tsv_to_ht.log', **add_args)
    ht = hl.import_table(args.input_file, no_header=True, impute=False, 
                         types={'f0':hl.tstr, 
                                'f1':hl.tint, # casting to non-str type doesn't work if there are empty strings (as there may be in some of the existing clump_results.txt files)
                                'f2':hl.tstr,
                                'f3':hl.tint, 
                                'f4':hl.tfloat,
                                'f5':hl.tint,
                                'f6':hl.tint,
                                'f7':hl.tint,
                                'f8':hl.tint,
                                'f9':hl.tint,
                                'f10':hl.tint,
                                'f11':hl.tstr})
    ht = ht.rename({'f0':'contig',
                    'f1':'F',
                    'f2':'varid',
                    'f3':'pos',
                    'f4':'P',
                    'f5':'TOTAL',
                    'f6':'NSIG',
                    'f7':'S05',
                    'f8':'S01',
                    'f9':'S001',
                    'f10':'S0001',
                    'f11':'SP2'})

    ht = ht.filter(ht.contig!='') # added because there were some files with extra lines with empty strings
    ht = ht.order_by('P')
    ht = ht.add_index()
    ht = ht.key_by(locus = hl.locus(contig=ht.contig, 
                                    pos=ht.pos,
                                    reference_genome='GRCh37'),
                   alleles = hl.array([ht.varid.split(':')[2],
                                       ht.varid.split(':')[3]])
    )
    ht = ht.annotate_globals(trait_type = args.trait_type,
                             phenocode = args.phenocode,
                             pheno_sex = args.pheno_sex,
                             coding = args.coding,
                             modifier = args.modifier)
    ht = ht.drop('contig','varid','pos')
    ht.describe()
#    ht.select().show()
    ht.write(args.output_file, overwrite=args.overwrite)

        
def export_ma_format(batch_size=256):
    r'''
    Export columns for .ma format (A1, A2, freq, beta, se, N) for select phenotypes
    '''
    meta_mt0 = hl.read_matrix_table('gs://ukb-diverse-pops/combined_results/meta_analysis.mt')
    
    highprev = hl.import_table(f'{ldprune_dir}/joined_ukbb_lancet_age_high_prev.tsv', impute=True)
    highprev = highprev.annotate(pheno = highprev.code.replace('_irnt',''))
    pheno_list = highprev.pheno.collect()
    pheno_list = [p for p in pheno_list if p is not None]
    meta_mt0 = meta_mt0.filter_cols(hl.literal(pheno_list).contains(meta_mt0.pheno))

    meta_mt0 = meta_mt0.annotate_cols(pheno_id = (meta_mt0.trait_type+'-'+
                                      meta_mt0.phenocode+'-'+
                                      meta_mt0.pheno_sex+
                                      hl.if_else(hl.len(meta_mt0.coding)>0, '-'+meta_mt0.coding, '')+
                                      hl.if_else(hl.len(meta_mt0.modifier)>0, '-'+meta_mt0.modifier, '')
                                      ).replace(' ','_').replace('/','_'))
    
    meta_mt0 = meta_mt0.annotate_rows(SNP = meta_mt0.locus.contig+':'+hl.str(meta_mt0.locus.position)+':'+meta_mt0.alleles[0]+':'+meta_mt0.alleles[1],
                                      A1 = meta_mt0.alleles[1], # .ma format requires A1 = effect allele, which in this case is A2 for UKB GWAS
                                      A2 = meta_mt0.alleles[0])

    meta_field_rename_dict = {'BETA':'b',
                          'SE':'se',
                          'Pvalue':'p',
                          'AF_Allele2':'freq',
                          'N':'N'}
    
    all_pops = ['AFR', 'AMR', 'CSA', 'EAS', 'EUR', 'MID']

    for pop in ['AFR','EUR']: #['AFR','AMR','CSA','EAS','EUR','MID']:
        print(f'not_{pop}')

        req_pop_list = [p for p in all_pops if p is not pop]

        loo_pop = meta_mt0.annotate_cols(idx = meta_mt0.meta_analysis_data.pop.index(hl.literal(req_pop_list))) # get index of which meta-analysis is the leave-on-out for current pop
        loo_pop = loo_pop.filter_cols(hl.is_defined(loo_pop.idx))
        
        annotate_dict = {}
        for field in ['AF_Allele2','BETA','SE','Pvalue','N']:
            annotate_dict.update({meta_field_rename_dict[field]: loo_pop.meta_analysis[field][loo_pop.idx]}) 
        batch_idx = 1

    export_out = f'{ldprune_dir}/loo/not_{pop}/batch{batch_idx}'
    while hl.hadoop_is_dir(export_out):
        batch_idx += 1
        export_out = f'{ldprune_dir}/loo/not_{pop}/batch{batch_idx}'
    checkpoint_path = f'gs://ukbb-diverse-temp-30day/loo/not_{pop}/batch{batch_idx}.mt'
#        print(f'\nCheckpointing to: {checkpoint_path}\n')
    loo_pop = loo_pop.checkpoint(checkpoint_path,
                                 _read_if_exists=True,
                                 overwrite=True)
    loo_pop = loo_pop.filter_entries(hl.is_defined(loo_pop.b))
    print(f'\nExporting to: {export_out}\n')
    hl.experimental.export_entries_by_col(mt = loo_pop,
                                          path = export_out,
                                          bgzip = True,
                                          batch_size = batch_size,
                                          use_string_key_as_file_name = True,
                                          header_json_in_file = False)
    

def mwzj_hts_by_tree(all_hts, temp_dir, globals_for_col_key, 
                     debug=False, inner_mode = 'overwrite', repartition_final: int = None,
                     read_if_exists = False):
    r'''
    Adapted from ukb_common mwzj_hts_by_tree()
    Uses read_clump_ht() instead of read_table()
    '''
    chunk_size = int(len(all_hts) ** 0.5) + 1
    outer_hts = []
    
    if read_if_exists: print('\n\nWARNING: Intermediate tables will not be overwritten if they already exist\n\n')
    
    checkpoint_kwargs = {inner_mode: True if not read_if_exists else False,
                         '_read_if_exists': read_if_exists} #
    if repartition_final is not None:
        intervals = ukb_common.get_n_even_intervals(repartition_final)
        checkpoint_kwargs['_intervals'] = intervals
    
    if debug: print(f'Running chunk size {chunk_size}...')
    for i in range(chunk_size):
        if i * chunk_size >= len(all_hts): break
        hts = all_hts[i * chunk_size:(i + 1) * chunk_size]
        if debug: print(f'Going from {i * chunk_size} to {(i + 1) * chunk_size} ({len(hts)} HTs)...')
        try:
            if isinstance(hts[0], str):
                def read_clump_ht(f):
                    ht = hl.read_table(f)
                    ht = ht.drop('idx')
#                    ht = ht.select(F=ht.F,
#                                   P=ht.P,
#                                   TOTAL=ht.TOTAL,
#                                   NSIG=ht.NSIG,
#                                   S05=ht.S05,
#                                   S01=ht.S01,
#                                   S001=ht.S001,
#                                   S0001=ht.S0001,
#                                   SP2=ht.SP2)
                    return ht
                hts = list(map(read_clump_ht, hts))
            ht = hl.Table.multi_way_zip_join(hts, 'row_field_name', 'global_field_name')
        except:
            if debug:
                print(f'problem in range {i * chunk_size}-{i * chunk_size + chunk_size}')
                _ = [ht.describe() for ht in hts]
            raise
        outer_hts.append(ht.checkpoint(f'{temp_dir}/temp_output_{i}.ht', **checkpoint_kwargs))
    ht = hl.Table.multi_way_zip_join(outer_hts, 'row_field_name_outer', 'global_field_name_outer')
    ht = ht.transmute(inner_row=hl.flatmap(lambda i:
                                           hl.cond(hl.is_missing(ht.row_field_name_outer[i].row_field_name),
                                                   hl.range(0, hl.len(ht.global_field_name_outer[i].global_field_name))
                                                   .map(lambda _: hl.null(ht.row_field_name_outer[i].row_field_name.dtype.element_type)),
                                                   ht.row_field_name_outer[i].row_field_name),
                                           hl.range(hl.len(ht.global_field_name_outer))))
    ht = ht.transmute_globals(inner_global=hl.flatmap(lambda x: x.global_field_name, ht.global_field_name_outer))
    mt = ht._unlocalize_entries('inner_row', 'inner_global', globals_for_col_key)
    return mt

def resume_mwzj(temp_dir, globals_for_col_key):
    ls = hl.hadoop_ls(temp_dir)
    paths = [x['path'] for x in ls if 'temp_output' in x['path'] ]
    chunk_size = len(paths)
    outer_hts = []
    for i in range(chunk_size):
        outer_hts.append(hl.read_table(f'{temp_dir}/temp_output_{i}.ht'))
    ht = hl.Table.multi_way_zip_join(outer_hts, 'row_field_name_outer', 'global_field_name_outer')
    ht = ht.transmute(inner_row=hl.flatmap(lambda i:
                                           hl.cond(hl.is_missing(ht.row_field_name_outer[i].row_field_name),
                                                   hl.range(0, hl.len(ht.global_field_name_outer[i].global_field_name))
                                                   .map(lambda _: hl.null(ht.row_field_name_outer[i].row_field_name.dtype.element_type)),
                                                   ht.row_field_name_outer[i].row_field_name),
                                           hl.range(hl.len(ht.global_field_name_outer))))
    ht = ht.transmute_globals(inner_global=hl.flatmap(lambda x: x.global_field_name, ht.global_field_name_outer))
    mt = ht._unlocalize_entries('inner_row', 'inner_global', globals_for_col_key)
    return mt

def join_clump_hts(pop, not_pop, high_quality=False, overwrite=False):
    r'''
    Wrapper for mwzj
    '''
    pop = pop.upper()
    pheno_manifest = hl.import_table(
    #            get_pheno_manifest_path(), 
            'gs://ukb-diverse-pops/ld_prune/phenotype_manifest.tsv.bgz', # hardcoded path to avoid having to change user-pays
            impute=True, 
            key=ukb_common.PHENO_KEY_FIELDS
            )
    pheno_manifest = pheno_manifest.annotate(pheno_id = pheno_manifest.filename.replace('.tsv.bgz',''))
    
    clump_results_dir = f'{ldprune_dir}/results{"_high_quality" if high_quality else ""}/{"not_" if not_pop else ""}{pop}'
    ls = hl.hadoop_ls(f'{clump_results_dir}/*')
    all_hts = [x['path'] for x in ls if 'clump_results.ht' in x['path']]
    
    temp_dir = f'gs://ukbb-diverse-temp-30day/nb-temp/{"not_" if not_pop else ""}{pop}{"-hq" if high_quality else ""}'
    globals_for_col_key = ukb_common.PHENO_KEY_FIELDS
    mt = mwzj_hts_by_tree(all_hts=all_hts,
                         temp_dir=temp_dir,
                         globals_for_col_key=globals_for_col_key)
#    mt = resume_mwzj(temp_dir=temp_dir, # NOTE: only use if all the temp hts have been created
#                     globals_for_col_key=globals_for_col_key)

    mt.write(get_clumping_results_path(pop=pop,
                                       not_pop=not_pop,
                                       high_quality=high_quality), 
             overwrite=overwrite)

def munge_mt(pop, not_pop, high_quality, parts=None):
    r'''
    For processing MTs before joining into the full clump mt
    '''
    pop_array = [p for p in all_pops if p!=pop] if not_pop else [pop]
    mt = hl.read_matrix_table(path=get_clumping_results_path(pop=pop,not_pop=not_pop,high_quality=high_quality),
                              _intervals=parts) # alternative path for results that are filtered to high quality variants before clumping
    print(f'{"not_" if not_pop else ""}{pop}: {mt.count_cols()}')
    mt = mt.annotate_cols(clump_pops = hl.literal(pop_array))
    if 'idx' in mt.entry.keys():
        mt = mt.select_entries(plink_clump=mt.entry.drop('idx')) # needed to clean up not_EAS and not_MID MTs (8/4/20)
    else:
        mt = mt.select_entries(plink_clump=mt.entry)
    mt.describe()
    return mt

def make_single_pop_clump_mt(high_quality=False, overwrite=False):
    mts = []
    not_pop=False
    
    n_partitions = 5000
    parts = hl.read_matrix_table(get_clumping_results_path(pop='EUR',not_pop=not_pop,high_quality=high_quality))._calculate_new_partitions(n_partitions)
    
    # for all single-pop clumping results
    for pop in all_pops:
        mts.append(munge_mt(pop=pop,
                            not_pop=not_pop,
                            high_quality=high_quality,
                            parts=parts))
    
    singlepop_mt = mts[0]
    for mt in mts[1:]:
        singlepop_mt = singlepop_mt.union_cols(mt, row_join_type='outer')
        
    singlepop_mt = singlepop_mt.collect_cols_by_key()
        
    singlepop_mt.write(get_clumping_results_path(pop='single_pop',high_quality=high_quality), 
                       overwrite=overwrite)
    # 2020-08-06 17:51:23 Hail: INFO: wrote matrix table with 24870911 rows and 7221 columns in 5000 partitions
    # Time: 2hr 48 min (ran with 2 workers, 150 preemptibles for first ~1hr, then removed all preemptibles)
    
def make_full_clump_mt(high_quality=False, overwrite=False):
    mts = []
    
    n_partitions = 5000 
    parts = hl.read_matrix_table(get_clumping_results_path(pop='EUR',not_pop=False,high_quality=high_quality))._calculate_new_partitions(n_partitions)

    for not_pop in [True, False]:
        for pop in all_pops:
            mts.append(munge_mt(pop=pop,
                                not_pop=not_pop,
                                high_quality=high_quality,
                                parts=parts))
    full_mt = mts[0]
    for mt in mts[1:]:
        full_mt = full_mt.union_cols(mt, row_join_type='outer')
        
    full_mt = full_mt.collect_cols_by_key()
    
    full_mt.write(get_clumping_results_path(pop='full',high_quality=high_quality), 
                  overwrite=overwrite)
    # 2020-08-07 06:17:10 Hail: INFO: wrote matrix table with 28024109 rows and 7221 columns in 5000 partitions
    # Time: 5hr 6 min (ran with 2 workers, 150 preemptibles for first ~1hr, then removed all preemptibles)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--n_threads', help='number of threads')
    parser.add_argument('--input_file', help='Input file of variant results')
    parser.add_argument('--output_file', help='Output file of variant results')
    parser.add_argument('--pop', type=str, help='Population to be left out')
    parser.add_argument('--trait_type', type=str, help='trait_type in meta-analyzed sumstats')
    parser.add_argument('--phenocode', type=str, help='phenocode in meta-analyzed sumstats')
    parser.add_argument('--pheno_sex', type=str, help='pheno_sex in meta-analyzed sumstats')
    parser.add_argument('--coding', type=str, default='', help='coding in meta-analyzed sumstats')
    parser.add_argument('--modifier', type=str, default='', help='modifier in meta-analyzed sumstats')
    parser.add_argument('--ht_to_tsv', action='store_true')
    parser.add_argument('--tsv_to_ht', action='store_true')
    parser.add_argument('--not_pop', action='store_true', help='whether pop set is a not_{pop}')
    parser.add_argument('--join_clump_hts', default=False, action='store_true')    
    parser.add_argument('--make_full_clump_mt', action='store_true')    
    parser.add_argument('--make_single_pop_clump_mt', action='store_true')    
    parser.add_argument('--batch_size', type=int, default=256, help='max number of phenotypes per batch for export_entries_by_col')
    parser.add_argument('--high_quality', default=False, action='store_true', help='Use high quality variants only')
    parser.add_argument('--overwrite', default=False, action='store_true', help='overwrite existing files')
    args = parser.parse_args()
    
#    try:        
    if args.ht_to_tsv:
        ht_to_tsv(args)
    elif args.tsv_to_ht:
        tsv_to_ht(args)
    elif args.join_clump_hts:
        join_clump_hts(pop=args.pop, 
                       not_pop=args.not_pop,
                       high_quality=args.high_quality,
                       overwrite=args.overwrite)
    elif args.make_full_clump_mt:
        make_full_clump_mt(high_quality=args.high_quality,
                           overwrite=args.overwrite)
    elif args.make_single_pop_clump_mt:
        make_single_pop_clump_mt(high_quality=args.high_quality, 
                                 overwrite=args.overwrite)

#    except:
#        hl.copy_log('gs://ukbb-diverse-temp-30day/nb_hail.log')
        