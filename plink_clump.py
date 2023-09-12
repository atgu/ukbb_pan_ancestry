#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 09:35:41 2020

Clumping GWAS results with PLINK

@author: nbaya
"""
from tkinter import FALSE
import hail as hl
import hailtop.batch as hb
from ukbb_pan_ancestry import bucket, POPS, PHENO_KEY_FIELDS
import argparse
import numpy as np
import pandas as pd

ldprune_dir = f'{bucket}/ld_prune'

all_chromosomes = list(range(1,23))+['X']

def read_plink_input_group_chrom(p, method, subset, chrom):
    r'''
    Reads input group of PLINK files into Batch.
    NOTE: This format allows for use of a single fam (not specific to a chromosome)
    '''
    assert method in {'clump','sbayesr'}
    
    if method == 'clump':
        subsets_dir = f'{ldprune_dir}/subsets_5k'
        prefix = f'{subsets_dir}/{subset}/{subset}'
        return p.read_input_group(bed=f'{prefix}.chr{chrom}.bed',
                                  bim=f'{prefix}.chr{chrom}.bim',
                                  fam=f'{prefix}.fam')
    elif method == 'sbayesr':
        subsets_dir = f'{ldprune_dir}/subsets_50k'
        prefix = f'{subsets_dir}/{subset}/{subset}'
        return p.read_input_group(bed=f'{prefix}.hm3.chr{chrom}.bed',
                                  bim=f'{prefix}.hm3.chr{chrom}.bim',
                                  fam=f'{prefix}.fam')

def get_pheno_list(pheno_manifest, pop: str, not_pop: bool, max_pops: bool, hq_phenos: bool):
    r'''
    Returns list of phenotypes in Pandas DataFrame `pheno_manifest` which include
    `pop` in the pops field of the pheno_manifest. 
    If `not_pop`=True, this will include all 6-population phenotypes and 
    phenotypes with num_pops=5 which do not have `pop` in the pops field.
    If `max_pops`=True, all phenotypes in `pheno_manifest` are returned. 
    ''' 
    assert not (not_pop and max_pops), '`not_pop` and `max_pops` cannot both be True'
    assert not (max_pops and hq_phenos), '`max_pops` and `hq_phenos` also cannot both be True'
    if not_pop and not hq_phenos:
        pheno_manifest = pheno_manifest[(pheno_manifest.num_pops==6)&
                                        ((pheno_manifest.num_pops==5)&(~pheno_manifest.pops.str.contains(pop)))]
    elif not max_pops and not hq_phenos: # if not not_pop and not full_pop
        pheno_manifest = pheno_manifest[pheno_manifest.pops.str.contains(pop)]
    elif not_pop and hq_phenos:
        pheno_manifest = pheno_manifest[pheno_manifest.pops_pass_qc.str.contains(pop, na=False)]
        pheno_manifest = pheno_manifest[(pheno_manifest.num_pops_pass_qc==3)|(pheno_manifest.num_pops_pass_qc==4)|(pheno_manifest.num_pops_pass_qc==5)]

    # pheno_manifest = pheno_manifest.drop_duplicates(subset='num_pops') # For testing purposes (to test single pheno for each unique num_pops)
    # pheno_manifest = pheno_manifest.drop_duplicates(subset='pops') # For testing purposes (to test single pheno for each unique pops set)
    
    # set fields according to arguments, makes it simpler when making pheno_list
    pheno_manifest['pop'] = pop
    pheno_manifest['not_pop'] = not_pop
    pheno_manifest['max_pops'] = max_pops
    pheno_manifest['hq_phenos'] = hq_phenos
    
    # .values --> Return a Numpy representation of the DataFrame
    pheno_manifest = pheno_manifest[pheno_manifest.num_pops_pass_qc > 0]

    if not max_pops and not hq_phenos:
        pheno_list = pheno_manifest[list(PHENO_KEY_FIELDS)+['pops','pop','not_pop','max_pops', 'hq_phenos']].values
    else:
        pheno_list = pheno_manifest[list(PHENO_KEY_FIELDS)+['pops_pass_qc','pop','not_pop','max_pops', 'hq_phenos']].values

#    seed = POPS.index(pop)#1
#    random.seed(a=seed)
#    random.shuffle(pheno_list)
#    
#    n_traits = 1
#    print(f'\nWARNING: For testing purposes, only using a max of {n_traits} traits (random seed: {seed})\n')
#    pheno_list = pheno_list[:n_traits]
    
#    print(pheno_list)
    print('\nNumber of phenotypes for '+
          ('max_pops=True' if max_pops else f'{"not_" if not_pop else ""}{pop}')+
          f': {len(pheno_list)}')

    return pheno_list

def get_pheno_id(trait_type, phenocode, pheno_sex, coding, modifier):
    r'''
    Creates unique ID for a phenotype
    '''
    pheno_id = (trait_type+'-'+
                str(phenocode)+'-'+ # TODO: remove str() for full pheno manifest
                pheno_sex+
                ('-'+coding if pd.notnull(coding) and len(coding)>0 else '')+
                ('-'+modifier if pd.notnull(modifier) and len(modifier)>0 else '')
                ).replace(' ','_').replace('/','_')
    return pheno_id

def get_pheno_key_dict(trait_type, phenocode, pheno_sex, coding, modifier, pops):
    r'''
    Used to make it simpler to pass phenotype information between methods.
    '''
    pheno_key_dict = {'trait_type': trait_type,
                      'phenocode': phenocode,
                      'pheno_sex': pheno_sex,
                      'coding': coding,
                      'modifier': modifier,
                      'pops': pops.split(',') if isinstance(pops, str) else pops}
    return pheno_key_dict

def get_sumstats(p, pop: str, not_pop: bool, max_pops: bool, hq_phenos: bool, pops: list, 
                 high_quality: bool, pheno_id: str, method: str, 
                 chromosomes: list = all_chromosomes):
    r'''
    Returns a dict of per-chromosome summary statistics output files.
    '''
    assert not (not_pop and max_pops), '`not_pop` and `max_pops` cannot both be True'
    assert not (max_pops and hq_phenos), '`max_pops` and `hq_phenos` also cannot both be True'
    assert method in {'clump','sbayesr'}
    if max_pops and len(pops)==1 and pop is None:
        pop = pops[0] # need to set this variable in order to find column indices later
        # pops is pheno_key_dict['pops']
    num_pops = len(pops)
    filename=f'{pheno_id}.tsv.bgz'
    trait_type = pheno_id.split('-')[0]
    trait_category = 'quant' if trait_type in ['continuous','biomarkers'] else 'binary'
    
    # change to variant manifest in public bucket 
    variant_manifest = p.read_input('gs://ukb-diverse-pops-public-free/sumstats_qc_analysis/full_variant_qc_metrics.txt.bgz')
    variant_manifest_tabix = p.read_input('gs://ukb-diverse-pops-public-free/sumstats_qc_analysis/full_variant_qc_metrics.txt.bgz.tbi')

    loo_6pop_dir = f'{ldprune_dir}/loo/sumstats/batch2'
    loo_6pop_ss_fname = f'{loo_6pop_dir}/{filename}'
    loo_6pop_tabix_fname = f'{loo_6pop_dir}_tabix/{filename}.tbi' # gs://ukb-diverse-pops/ld_prune/loo/sumstats/batch2_tabix

    loo_3pop_dir = f'{ldprune_dir}/loo/sumstats/hq/221215_hq_loo_3pop' # gs://ukb-diverse-pops/ld_prune/
    loo_3pop_ss_fname = f'{loo_3pop_dir}/{filename}'
    loo_3pop_tabix_fname = f'{loo_3pop_dir}_tabix/{filename}.tbi'

    ss_dir = f'{bucket}/sumstats_flat_files'
    ss_fname = f'{ss_dir}/{filename}'
    tabix_fname = f'{ss_dir}_tabix/{filename}.tbi'
    
    get_ss = p.new_job(name=f'get_ss_{pheno_id}')
    get_ss = get_ss.image('gcr.io/ukbb-diversepops-neale/nbaya_tabix:latest')
    get_ss.storage('100M') # default: 1G
    get_ss.cpu(1)
    bgz_fname = f'{get_ss.ofile}.bgz'
    tbi_fname = f'{get_ss.ofile}.bgz.tbi'
    get_ss.command('set -ex')
    variant_manifest_bgz = f'{get_ss.ofile}.variants.bgz'
    variant_manifest_tbi = f'{get_ss.ofile}.variants.bgz.tbi'
    get_ss.command(' '.join(['mv',variant_manifest, variant_manifest_bgz]))
    get_ss.command(' '.join(['mv',variant_manifest_tabix, variant_manifest_tbi]))
    
    if not_pop and not hq_phenos and hl.hadoop_is_file(loo_6pop_ss_fname) and hl.hadoop_is_file(loo_6pop_tabix_fname): # phenotype is 6-pop and has leave-one-out sumstats generated. 
#        assert False, "don't run 6-pop LOO"
        print(f'Using 6-pop LOO sumstats for {pheno_id} ({"not_" if not_pop else ""}{pop})')
        ss = p.read_input(loo_6pop_ss_fname)
        tabix = p.read_input(loo_6pop_tabix_fname)

        get_ss.command(' '.join(['mv',ss,bgz_fname])) # necessary instead of changing path extension for input files
        get_ss.command(' '.join(['mv',tabix,tbi_fname])) # necessary instead of changing path extension for input files
        get_ss.command('\n'.join(
                f'''
                tabix -h {bgz_fname} {chrom} | \
                cut -f5,{6+POPS.index(pop)} | \
                sed 's/pval_not_{pop}/P/g' | \
                awk '$2!="NA" {{print}}' > {get_ss[f'ofile_{chrom}']}
                '''
                for chrom in chromosomes
                )
        )

    elif not_pop and hq_phenos and hl.hadoop_is_file(loo_3pop_ss_fname) and hl.hadoop_is_file(loo_3pop_tabix_fname):
        print(f'Using {num_pops}-pop LOO sumstats for {pheno_id} (not_{pop})')
        ss = p.read_input(loo_3pop_ss_fname)
        tabix = p.read_input(loo_3pop_tabix_fname)

        get_ss.command(' '.join(['mv',ss,bgz_fname])) # necessary instead of changing path extension for input files
        get_ss.command(' '.join(['mv',tabix,tbi_fname])) # necessary instead of changing path extension for input files
        
        awk_arg2 = '$2!="NA"' + (' && $3!="false"' if high_quality else '')

        if high_quality:
            get_ss.command('\n'.join(
                f'''paste <( tabix {bgz_fname} {chrom} | \
                    cut -f3,{4+pops.index(pop)}) \
                    <( tabix {variant_manifest_bgz} {chrom} | \
                    awk '{{ print $9 }}' ) | \
                    awk -v OFS="\t" '{{if({awk_arg2}) print $1,10^-$2}}' | awk 'NR == 1 {{print "SNP\tP"}}1' > {get_ss[f"ofile_{chrom}"]}
                    '''
                    for chrom in chromosomes
                    )
        )

        else:
            get_ss.command('\n'.join(f'''
                    tabix {bgz_fname} {chrom} | \
                    cut -f3,{4+pops.index(pop)} | \
                        awk -F'\t' 'OFS="\t" {{print $1,10^-$2}}' | \
                        awk 'NR == 1 {{print "SNP\tP"}}1' | \
                            awk '$2!="NA" {{print}}' > {get_ss[f'ofile_{chrom}']}
                    '''
                    for chrom in chromosomes
                    )
            )       

    elif hl.hadoop_is_file(ss_fname) and hl.hadoop_is_file(tabix_fname): # this conditional block must come after checking for 6-pop LOO results
        print(f'Using {num_pops}-pop sumstats for {pheno_id} '+
              (f'({"not_" if not_pop else ""}{pop})' if not max_pops else '(max_pops=True)'))
        ss = p.read_input(ss_fname)
        tabix = p.read_input(tabix_fname)
    
        get_ss.command(' '.join(['mv',ss,bgz_fname])) # necessary instead of changing path extension for input files
        get_ss.command(' '.join(['mv',tabix,tbi_fname])) # necessary instead of changing path extension for input files
        
        if not_pop or (max_pops and len(pops)>1):
            pval_col_idx = 8 if trait_category == 'quant' else 9 # due to additional AF columns in binary traits, pvalue column location may change
            awk_arg1 = ''
            awk_arg2 = '$2!="NA"'+ ('&& $3!="false"' if high_quality else '') # exclude pval(col 2)=NA; if high_quality: exclude high_quality(col 3)=false
        else: # if clumping single population results
            pval_col_idx = ((9+(trait_category == 'binary'))+ 
                            ((4+(trait_category =='binary')+1) if num_pops>1 else 0)+
                            ((trait_category == 'binary')+3)*num_pops+
                             pops.index(pop)+1)

            low_confidence_col_idx = ((9+(trait_category == 'binary'))+ # first 4 cols
                                      ((4+(trait_category =='binary')+1) if num_pops>1 else 0)+ # meta-analysis fields
                                      ((trait_category == 'binary')+4)*num_pops+ # per-pop fields
                                      pops.index(pop)+1)
            awk_arg1 = f', $3=${low_confidence_col_idx}'
            awk_arg2 = '$2!="NA" && $3!="true"' + (' && $4!="false"' if high_quality else '') # exclude pval(col 2)=NA, low_confidence(col 3)=True; if high_quality: exclude high_quality(col 4)=False
        
        # TODO: If possible, consolidate the following blocks
        if high_quality: 
            get_ss.command('\n'.join(
                    f'''
                    paste <( tabix {bgz_fname} {chrom} | \
                            awk '{{print $1=$1":"$2":"$3":"$4, $2=${pval_col_idx}{awk_arg1}}}') \
                          <( tabix {variant_manifest_bgz} {chrom} | \
                            awk '{{ $9 }}' ) | \
                    awk -v OFS="\t" '{{if({awk_arg2}) print $1,10^-$2}}' | awk 'NR == 1 {{print "SNP\tP"}}1' > {get_ss[f"ofile_{chrom}"]}
                    '''
                    for chrom in chromosomes
                    )
            )
        else: 
            get_ss.command('\n'.join(
                    f'''
                    tabix -h {bgz_fname} {chrom} | \
                    awk '{{print $1=$1":"$2":"$3":"$4, $2=${pval_col_idx}{awk_arg1}}}') | \
                    awk '{{if({awk_arg2}) print $1,$1,10^-$2}}' | awk 'NR == 1 {{print "SNP\tP"}}1' > {get_ss[f"ofile_{chrom}"]}
                    '''
                    for chrom in chromosomes
                    )
            )
                              

    ss_dict = {
            chrom: get_ss[f'ofile_{chrom}']
            for chrom in chromosomes
    }
        
    # if filter_hm3:
    #     filter = p.new_python_job(name='Filtering to HM3 SNPs')
    #     ss_dict_filtered = filter.call(get_hapmap3_snps, hm3_snp_list = '', ss = ss_dict) ## can this be get_ss or does it have to be a new_python_job?? like "j = format_b.new_python_job(name=f'Formatting: {sst}')"" (but format_b would be p in this case) then "j.call(format_input, sst, SNP_name, A1_name, A2_name, beta, pval, pheno, args.out_dir, args.snp_info)"
    #     return ss_dict_filtered
    
    return ss_dict
    
def get_adj_betas(p, pop, not_pop, max_pops, hq_phenos, pheno_key_dict, pheno_id, high_quality, 
                  hail_script,clump_hm3: bool=False):
    r'''
    Wrapper method for both PLINK clumping and SBayesR
    '''

    output_dir = (f'{ldprune_dir}/results{"_hq_phenos" if hq_phenos else ""}'+('_22115/' if not clump_hm3 else '_HM3/')+
                ('max_pops' if max_pops else f'{"not_" if not_pop else ""}{pop}')+('_hq' if high_quality else "")+
                f'/{pheno_id}')


    clump_output_txt = f'{output_dir}/clump_results.txt' # PLINK clump output txt file
    clump_output_ht = f'{output_dir}/clump_results.ht' # PLINK clump output hail table
    
#    sbayesr_output_txt = f'{output_dir}/sbayesr_results-test.txt' # SBayesR output txt file
#    sbayesr_output_ht = f'{output_dir}/sbayesr_results-test.ht' # SBayesR output hail table
    
    # TODO: make argument
    overwrite = True
    
    clump_file_exists = hl.hadoop_is_file(f'{clump_output_ht}/_SUCCESS')
    if not clump_file_exists or overwrite:
        if clump_file_exists and overwrite:
            print(f'\n\nWARNING: Existing results will be overwritten for {pheno_id} in {output_dir}!\n')

        ss_dict = get_sumstats(p=p,
                            pop=pop, 
                            not_pop=not_pop,
                            max_pops=max_pops,
                            hq_phenos=hq_phenos,
                            pops=pheno_key_dict['pops'],
                            high_quality=high_quality,
                            pheno_id=pheno_id,
                            method='clump')

        run_method(p=p, 
                   pop=pop, 
                   not_pop=not_pop,
                   max_pops=max_pops,
                   pheno_id=pheno_id, 
                   pheno_key_dict=pheno_key_dict,
                   hail_script=hail_script, 
                   output_txt=clump_output_txt, 
                   output_ht=clump_output_ht,
                   ss_dict=ss_dict,
                   method='clump',
                   hm3=clump_hm3)
    else:
        print(f'\n\nSkipping {pheno_id} because results ht exists and overwrite=False\n')
        
#    if not hl.hadoop_is_file(f'{sbayesr_output_ht}/_SUCCESS'):
#        if overwrite:
#            print('\n\nWARNING: Existing results will be overwritten!\n')
#            
#        meta_ss = p.read_input(f'{ldprune_dir}/loo/not_AFR/batch2/biomarkers-30740-30740.tsv.bgz')
#
#        run_method(p=p, 
#                   pop=pop, 
#                   pheno=pheno, 
#                   coding=coding, 
#                   trait_type=trait_type, 
#                   hail_script=hail_script, 
#                   output_txt=sbayesr_output_txt,
#                   output_ht=sbayesr_output_ht,
#                   ss=meta_ss,
#                   method='sbayesr')
        

def run_method(p, pop, not_pop, max_pops, pheno_key_dict, pheno_id, 
               hail_script, output_txt, output_ht, ss_dict, method, hm3: bool=False):
    r'''
    Runs either PLINK clump (method = 'clump') or SBayesR (method = 'sbayesr')
    '''
    assert method in {'clump','sbayesr'}
    
    task_suffix = (f'{"not_" if not_pop else ""}{pop}' if not max_pops else 'max_pops')+f'-{pheno_id}'
    # TODO: if method = 'sbayesr' check if LD matrix has already been calculated
    
    tasks = []
    
    ref_subset = '-'.join(pheno_key_dict['pops'] if max_pops else ([p for p in pheno_key_dict['pops'] if p != pop] if not_pop else [pop]))
    print(f'Using LD reference panel of {ref_subset}')
        
    ## run plink clumping
    for chrom, ss_chrom in ss_dict.items():
        
        ## read ref ld plink files 
        bfile = read_plink_input_group_chrom(p=p, 
                                             method=method,
                                             subset=ref_subset,
                                             chrom=chrom)

        ## read HM3 SNP list
        if hm3:
            hm3_list = p.read_input(f'{bucket}/ktsuo_unrelateds_tmp/hm3_hg37.tsv')
        
        get_betas = p.new_job(name=f'{method}_{task_suffix}_chr{chrom}')
        
        # TODO: change image to include GCTB if running SBayesR?
        get_betas.cpu(1) # plink clump cannot multithread
        
        get_betas.command('set -ex')
        
        if method == 'clump' and hm3:
            file_info = hl.utils.hadoop_stat(f'{bucket}/ld_prune/subsets_5k/{ref_subset}/{ref_subset}.chr{chrom}.bed')
            size_bytes = file_info['size_bytes']
            size_gigs = size_bytes / (1024 * 1024 * 1024)
            job_storage = round(10.0 + 2.0 * size_gigs)
            get_betas.storage(job_storage) # default: 5G
            # clump_memory = -15*(chrom-1)+400 # Memory requested for PLINK clumping in MB. equation: -15*(chrom-1) + 500 is based on 400 MB for chr 1, 80 MB for chr 22
            clump_memory = 4 # in GB
            get_betas.memory(f'{clump_memory}Gi') # default: 30G
            get_betas.command(f'head {ss_chrom}')
            get_betas.command(' '.join(['plink',
                                        '--bfile', str(bfile),
                                        '--memory',str(clump_memory*1000), # memory in MB
                                        '--threads','1', # explicitly set threads to 1
                                        '--clump', ss_chrom,
                                        '--extract', str(hm3_list),
                                        '--clump-field P',
                                        '--clump-snp-field SNP',
                                        '--clump-p1 1',
                                        '--clump-p2 1',
                                        '--clump-r2 0.1',
                                        '--clump-kb 500',
                                        '--maf 0.01',
                                        '--geno 0.05',
                                        '--hwe 1e-6',
                                        '--output-chr M', # necessary to code chr X as 'X' instead of '23', which isn't allowed as a contig in Hail's GRCh37 locus
                                        '--chr', str(chrom),
                                        '--out',f'{get_betas.ofile}_tmp']))
            get_betas.command(' '.join(['awk',"'{ print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12 }'",
                                        "OFS='\t'",f'{get_betas.ofile}_tmp.clumped',
                                        '2>','/dev/null','|',
                                        'tail -n+2','|',  # don't include header
                                        "sed '/^[[:space:]]*$/d'", # remove 2 empty lines created by PLINK at the end of the output file 
                                        '>',str(get_betas.ofile)]))
        elif method == 'clump' and not hm3:
            file_info = hl.utils.hadoop_stat(f'{bucket}/ld_prune/subsets_5k/{ref_subset}/{ref_subset}.chr{chrom}.bed')
            size_bytes = file_info['size_bytes']
            size_gigs = size_bytes / (1024 * 1024 * 1024)
            job_storage = round(10.0 + 2.0 * size_gigs)
            get_betas.storage(job_storage) # default: 5G
            # clump_memory = -15*(chrom-1)+400 # Memory requested for PLINK clumping in MB. equation: -15*(chrom-1) + 500 is based on 400 MB for chr 1, 80 MB for chr 22
            clump_memory = 4 # in GB
            get_betas.memory(f'{clump_memory}Gi') # default: 30G
            get_betas.command(f'head {ss_chrom}')
            get_betas.command(' '.join(['plink',
                                        '--bfile', str(bfile),
                                        '--memory',str(clump_memory*1000), # memory in MB
                                        '--threads','1', # explicitly set threads to 1
                                        '--clump-field P',
                                        '--clump-snp-field SNP',
                                        '--clump-p1 1',
                                        '--clump-p2 1',
                                        '--clump-r2 0.1',
                                        '--clump-kb 500',
                                        '--output-chr M', # necessary to code chr X as 'X' instead of '23', which isn't allowed as a contig in Hail's GRCh37 locus
                                        '--chr', str(chrom),
                                        '--out',f'{get_betas.ofile}_tmp']))
            get_betas.command(' '.join(['awk',"'{ print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12 }'",
                                        "OFS='\t'",f'{get_betas.ofile}_tmp.clumped',
                                        '2>','/dev/null','|',
                                        'tail -n+2','|',  # don't include header
                                        "sed '/^[[:space:]]*$/d'", # remove 2 empty lines created by PLINK at the end of the output file 
                                        '>',str(get_betas.ofile)]))
        elif method == 'sbayesr':
            ldm_type = 'full' # options: full, sparse
            ldm_path = f'{ldprune_dir}/subsets_50k/not_{pop}/ldm/not_{pop}.hm3.chr{chrom}.maf_gt_0.ldm.{ldm_type}'
#            ldm_path = f'{ldprune_dir}/subsets_50k/not_{pop}/ldm/not_{pop}.hm3.chr{chrom}.maf_gt_0.chisq_5.ldm.sparse'
            
            if hl.hadoop_is_file(f'{ldm_path}.info') and hl.hadoop_is_file(f'{ldm_path}.bin'):
                ldm = p.read_input_group(info=f'{ldm_path}.info',
                                         bin=f'{ldm_path}.bin')
            else:
                make_ldm = p.new_job(name=f'make_{ldm_type}_ldm_{task_suffix}.chr{chrom}')
                make_ldm.memory('60G')
                make_ldm.command(' '.join(['wget',
                                        'https://cnsgenomics.com/software/gctb/download/gctb_2.0_Linux.zip',
                                        '-P', '~/']))
                make_ldm.command(' '.join(['unzip','~/gctb_2.0_Linux.zip','-d','~/']))
                make_ldm.command(' '.join(['ls','-ltrR','~/']))
                make_ldm.command(' '.join(['mv','~/gctb_2.0_Linux/gctb','/usr/local/bin/']))
                make_ldm.command(' '.join(['plink',
                                           '--bfile', str(bfile),
                                           '--maf 0.0000000001',
                                           '--make-bed',
                                           '--out', f'{make_ldm.ofile}_tmp1']))
                make_ldm.command(' '.join(['gctb',
                                            '--bfile', f'{make_ldm.ofile}_tmp1',
#                                            '--snp 1-1000',
                                            f'--make-{ldm_type}-ldm', 
                                            '--out',f'{make_ldm.ofile}_tmp2']))
                # TODO: use both .bin and .info files
                make_ldm.command(' '.join(['mv',f'{make_ldm.ofile}_tmp2.ldm.{ldm_type}', str(make_ldm.ofile)]))
                p.write_output(make_ldm.ofile, ldm_path)
                ldm = make_ldm.ofile
                
            get_betas.declare_resource_group(out={'log':'{root}.log',
                                                  'snpRes':'{root}.snpRes',
                                                  'parRes':'{root}.parRes',
                                                  'mcmcsamples.SnpEffects':'{root}.mcmcsamples.SnpEffects',
                                                  'mcmcsamples.Par':'{root}.mcmcsamples.Par'})
            get_betas.command(' '.join(['wget',
                                        'https://cnsgenomics.com/software/gctb/download/gctb_2.0_Linux.zip',
                                        '-P', '~/']))
            get_betas.memory('18G')
            get_betas.command(' '.join(['unzip','~/gctb_2.0_Linux.zip','-d','~/']))
            get_betas.command(' '.join(['ls','-ltrR','~/']))
            get_betas.command(' '.join(['mv','~/gctb_2.0_Linux/gctb','/usr/local/bin/']))
            get_betas.command(' '.join(['gctb',
                                        '--sbayes R', 
                                        '--ldm', str(ldm),
                                        '--pi 0.95,0.02,0.02,0.01',
                                        '--gamma 0.0,0.01,0.1,1',
                                        '--gwas-summary', f' <( gunzip -c {ss_chrom} | grep -v "NA" )',
                                        '--chain-length 10000',
                                        '--burn-in 2000',
                                        '--out-freq 10',
                                        '--out',f'{get_betas.out}']))
            get_betas.command(' '.join(['head',f'{get_betas.out}.snpRes']))
            get_betas.command(' '.join(['mv',f'{get_betas.out}.snpRes', str(get_betas.ofile)]))
            
        tasks.append(get_betas)

    get_betas_sink = p.new_job(name=f'{method}_sink_{task_suffix}')
    get_betas_sink.command(f'cat {" ".join([t.ofile for t in tasks])} > {get_betas_sink.ofile}') # this task implicitly depends on the chromosome scatter tasks
    p.write_output(get_betas_sink.ofile, output_txt)
    
    ## import as hail table and save
    n_threads = 8
    tsv_to_ht = p.new_job(name=f'{method}_to_ht_{task_suffix}')
    tsv_to_ht = tsv_to_ht.image('gcr.io/ukbb-diversepops-neale/nbaya_hail:latest')
    tsv_to_ht.storage('1G')
    tsv_to_ht.memory('100M')
    tsv_to_ht.cpu(n_threads)
    tsv_to_ht.depends_on(get_betas_sink)
    tsv_to_ht.command('set -ex')
    tsv_to_ht.command(' '.join(['PYTHONPATH=$PYTHONPATH:/',
                        'PYSPARK_SUBMIT_ARGS="--conf spark.driver.memory=4g --conf spark.executor.memory=24g pyspark-shell"']))
    tsv_to_ht.command(' '.join(['python3', str(hail_script),
                          '--input_file', f'"{output_txt}"', # output_txt must be doubly enclosed by quotes needed for files with "|" in their pheno_id
                          '--tsv_to_ht',
                          '--trait_type', f'''"{pheno_key_dict['trait_type']}"''',
                          '--phenocode', f'''"{pheno_key_dict['phenocode']}"''',
                          '--pheno_sex', f'''"{pheno_key_dict['pheno_sex']}"''',
                          '--output_file', f'"{output_ht}"', # output_ht must be doubly enclosed by quotes needed for files with "|" in their pheno_id
                          '--overwrite']+
                          (['--coding',f'''"{pheno_key_dict['coding']}"'''] if pheno_key_dict['coding'] != '' else [])+
                          (['--modifier',f'''"{pheno_key_dict['modifier']}"'''] if pheno_key_dict['modifier'] != '' else [])))

def main(args):    
    backend = hb.ServiceBackend(billing_project='ukb_diverse_pops', bucket='ukbb-diverse-temp-30day/kt-batch-tmp')
#    backend = batch.LocalBackend(tmp_dir='/tmp/batch/')

    high_quality = args.high_quality
    not_pop = args.not_pop
    max_pops = args.max_pops
    hq_phenos = args.hq_phenos
    clump_hm3 = args.clump_hm3

    paridx = 0
    parsplit = 1
    
    p = hb.Batch(name=(f'clump{"-hq" if high_quality else ""}{"-max_pops" if max_pops else ""}{"-hq_phenos" if hq_phenos else ""}'+
                             f'-{paridx}-{parsplit}'), 
                       backend=backend, default_image='gcr.io/ukbb-diversepops-neale/nbaya_plink:0.1',
                       default_storage='500Mi', default_cpu=8)

    # download hail script to VM
    hail_script = p.read_input(f'{ldprune_dir}/scripts/python/plink_clump_hail.py')
    
    # get phenotype manifest   
    pheno_manifest = pd.read_table('gs://ukb-diverse-pops/combined_results/221215_phenotype_manifest.tsv.bgz', header=0, sep='\t', compression='gzip')
    pheno_manifest.fillna('', inplace=True)

    # max_pops returns ALL phenotypes, but many phenos don't have pops_pass_qc, so need to filter first to hq_phenos
    if max_pops:
        pheno_list = get_pheno_list(pheno_manifest=pheno_manifest,
                                    pop=None,
                                    not_pop=False,
                                    max_pops=True,
                                    hq_phenos=False)

    elif hq_phenos:
        pheno_list = []
        pops_noEUR = ['AFR', 'AMR', 'CSA', 'EAS', 'MID']

        for pop in pops_noEUR:
            pheno_list.append(get_pheno_list(pheno_manifest=pheno_manifest,
                                        pop=pop,
                                        not_pop=True, 
                                        max_pops=False,
                                        hq_phenos=True))

    elif not_pop:
        pheno_list = []
        # for not_pop in [True, False]: 
        for pop in POPS:
            pheno_list.append(get_pheno_list(pheno_manifest=pheno_manifest,
                                                pop=pop,
                                                not_pop=True,
                                                max_pops=False,
                                                hq_phenos=False))

    else:
        pheno_list = []
        for pop in POPS:
                pheno_list.append(get_pheno_list(pheno_manifest=pheno_manifest,
                                                 pop=pop,
                                                 not_pop=False,
                                                 max_pops=False,
                                                 hq_phenos=False))

    
    # if hq_phenos:
    if not max_pops: # need for hq_phenos and single pops
        pheno_list = np.concatenate(pheno_list) # necessary for hq_phenos because getting an array of arrays per not_POP
    idx_to_run = range(paridx, len(pheno_list), parsplit)
    pheno_list = [pheno_info for idx, pheno_info in enumerate(pheno_list) if idx in idx_to_run]
    print(f'Number of phenos: {len(pheno_list)}')
    
    if max_pops or hq_phenos or not_pop:
        for trait_type, phenocode, pheno_sex, coding, modifier, pops_pass_qc, pop, not_pop, max_pops, hq_phenos in pheno_list:
            
            pheno_key_dict = get_pheno_key_dict(trait_type, phenocode, pheno_sex, coding, modifier, pops_pass_qc)
            pheno_id = get_pheno_id(trait_type, phenocode, pheno_sex, coding, modifier)

            get_adj_betas(p=p,
                        pop=pop, 
                        not_pop=not_pop,
                        max_pops=max_pops,
                        hq_phenos=hq_phenos,
                        pheno_key_dict=pheno_key_dict,
                        high_quality=high_quality,
                        pheno_id=pheno_id, 
                        hail_script=hail_script,
                        clump_hm3=clump_hm3)
    else:
        for trait_type, phenocode, pheno_sex, coding, modifier, pops, pop, not_pop, max_pops, hq_phenos in pheno_list:
            
            pheno_key_dict = get_pheno_key_dict(trait_type, phenocode, pheno_sex, coding, modifier, pops)
            pheno_id = get_pheno_id(trait_type, phenocode, pheno_sex, coding, modifier)

            get_adj_betas(p=p,
                        pop=pop, 
                        not_pop=not_pop,
                        max_pops=max_pops,
                        hq_phenos=hq_phenos,
                        pheno_key_dict=pheno_key_dict,
                        high_quality=high_quality,
                        pheno_id=pheno_id, 
                        hail_script=hail_script,
                        clump_hm3=clump_hm3)

            
    # p.run(open=True)
    p.run()
    
#    if type(backend)==batch.ServiceBackend:
#        print('running')
#    else:
#        p.run(verbose=True,
#              delete_scratch_on_exit=True)
    
    # backend.close()


if __name__=="__main__":

    parser = argparse.ArgumentParser()

    # parser.add_argument('--not_pop', action='store_true', help='whether pop set is a not_pop')
    parser.add_argument('--max_pops', action='store_true', help='Include to get "max_pops" clumping results')
    parser.add_argument('--hq_phenos', action='store_true', help='Include to get clumping results for hq phenotypes only')
    parser.add_argument('--high_quality', action='store_true', help='Include to clump based on high quality variants only')
    parser.add_argument('--not_pop', action='store_true')
    parser.add_argument('--clump_hm3', action='store_true', help='Include to clump based on HM3 SNPs only')

    args = parser.parse_args()    


    # hl.init(default_reference='GRCh38', 
    #         spark_conf={'spark.hadoop.fs.gs.requester.pays.mode': 'AUTO', 
    #                     'spark.hadoop.fs.gs.requester.pays.project.id': 'ukbb-diversepops-neale'})

    main(args)
