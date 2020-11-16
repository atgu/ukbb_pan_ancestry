#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 09:35:41 2020

Clumping GWAS results with PLINK

@author: nbaya
"""

import hail as hl
import hailtop.batch as hb
from ukbb_pan_ancestry import bucket, POPS, PHENO_KEY_FIELDS

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

def get_pheno_list(pheno_manifest, pop: str, not_pop: bool, max_pops: bool):
    r'''
    Returns list of phenotypes in Pandas DataFrame `pheno_manifest` which include
    `pop` in the pops field of the pheno_manifest. 
    If `not_pop`=True, this will include all 6-population phenotypes and 
    phenotypes with num_pops=5 which do not have `pop` in the pops field.
    If `max_pops`=True, all phenotypes in `pheno_manifest` are returned.
    ''' 
    assert not (not_pop and max_pops), '`not_pop` and `max_pops` cannot both be True'
    if not_pop:
        pheno_manifest = pheno_manifest[(pheno_manifest.num_pops==6)&
                                        ((pheno_manifest.num_pops==5)&(~pheno_manifest.pops.str.contains(pop)))]
    elif not max_pops: # if not not_pop and not full_pop
        pheno_manifest = pheno_manifest[pheno_manifest.pops.str.contains(pop)]

    # pheno_manifest = pheno_manifest.drop_duplicates(subset='num_pops') # For testing purposes (to test single pheno for each unique num_pops)
    # pheno_manifest = pheno_manifest.drop_duplicates(subset='pops') # For testing purposes (to test single pheno for each unique pops set)
    
    # set fields according to arguments, makes it simpler when making pheno_list
    pheno_manifest['pop'] = pop
    pheno_manifest['not_pop'] = not_pop
    pheno_manifest['max_pops'] = max_pops
    
    pheno_list = pheno_manifest[list(PHENO_KEY_FIELDS)+['pops','pop','not_pop','max_pops']].values

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
                phenocode+'-'+
                pheno_sex+
                ('-'+coding if len(coding)>0 else '')+
                ('-'+modifier if len(modifier)>0 else '')
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

def get_sumstats(p, pop: str, not_pop: bool, max_pops: bool, pops: list, 
                 high_quality: bool, pheno_id: str, method: str, 
                 chromosomes: list = all_chromosomes):
    r'''
    Returns a dict of per-chromosome summary statistics output files.
    '''
    assert not (not_pop and max_pops), '`not_pop` and `max_pops` cannot both be True'
    assert method in {'clump','sbayesr'}
    if max_pops and len(pops)==1 and pop is None:
        pop = pops[0] # need to set this variable in order to find column indices later
    num_pops = len(pops)
    filename=f'{pheno_id}.tsv.bgz'
    trait_type = pheno_id.split('-')[0]
    trait_category = 'quant' if trait_type in ['continuous','biomarkers'] else 'binary'
    
    variant_manifest = p.read_input(f'{ldprune_dir}/variant_qc/full_variant_qc_metrics.txt.bgz')
    variant_manifest_tabix = p.read_input(f'{ldprune_dir}/variant_qc_tabix/full_variant_qc_metrics.txt.bgz.tbi')

    loo_6pop_dir = f'{ldprune_dir}/loo/sumstats/batch2'
    loo_6pop_ss_fname = f'{loo_6pop_dir}/{filename}'
    loo_6pop_tabix_fname = f'{loo_6pop_dir}_tabix/{filename}.tbi'
    
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
        
    if not_pop and hl.hadoop_is_file(loo_6pop_ss_fname) and hl.hadoop_is_file(loo_6pop_tabix_fname): # phenotype is 6-pop and has leave-one-out sumstats generated. 
#        assert False, "don't run 6-pop LOO"
        print(f'Using 6-pop LOO sumstats for {pheno_id} ({"not_" if not_pop else ""}{pop})')
        ss = p.read_input(loo_6pop_ss_fname)
        tabix = p.read_input(loo_6pop_tabix_fname)
        
        get_ss.command(' '.join(['mv',ss,bgz_fname])) # necessary instead of changing path extension for input files
        get_ss.command(' '.join(['mv',tabix,tbi_fname])) # necessary instead of changing path extension for input files
        get_ss.command('\n'.join(
                f'''
                tabix -h {bgz_fname} {chrom} | \\
                cut -f5,{6+POPS.index(pop)} | \\
                sed 's/pval_not_{pop}/P/g' | \\
                awk '$2!="NA" {{print}}' > {get_ss[f'ofile_{chrom}']}
                '''
                for chrom in chromosomes
                )
        )
    elif hl.hadoop_is_file(ss_fname) and hl.hadoop_is_file(tabix_fname): # this conditional block must come after checking for 6-pop LOO results
        print(f'Using {num_pops}-pop sumstats for {pheno_id} '+
              ('({"not_" if not_pop else ""}{pop})' if not max_pops else '(max_pops=True)'))
        ss = p.read_input(ss_fname)
        tabix = p.read_input(tabix_fname)
    
        get_ss.command(' '.join(['mv',ss,bgz_fname])) # necessary instead of changing path extension for input files
        get_ss.command(' '.join(['mv',tabix,tbi_fname])) # necessary instead of changing path extension for input files
        
        if not_pop or (max_pops and len(pops)>1): # if clumping on meta-analyzed sumstats
            pval_col_idx = 8 if trait_category == 'quant' else 9 # due to additional AF columns in binary traits, pvalue column location may change
            awk_arg1 = ''
            awk_arg2 = '$2!="NA"'+ ('&& $3!="false"' if high_quality else '') # exclude pval(col 2)=NA; if high_quality: exclude high_quality(col 3)=false
            sed_arg = "-e 's/pval_meta/P/g'"
        else: # if clumping single population results
            pval_col_idx = (4+((4+(trait_category =='binary')+1) if num_pops>1 else 0)+
                            ((trait_category == 'binary')+3)*num_pops+
                             pops.index(pop)+1)

            low_confidence_col_idx = (4+ # first 4 cols
                                      ((4+(trait_category =='binary')+1) if num_pops>1 else 0)+ # meta-analysis fields
                                      ((trait_category == 'binary')+4)*num_pops+ # per-pop fields
                                      pops.index(pop)+1)
            awk_arg1 = f', $3=${low_confidence_col_idx}'
            awk_arg2 = '$2!="NA" && $3!="true"' + (' && $4!="false"' if high_quality else '') # exclude pval(col 2)=NA, low_confidence(col 3)=True; if high_quality: exclude high_quality(col 4)=False
            sed_arg =  f"-e 's/pval_{pop}/P/g'" # sed argument for replacing column name
        
        # TODO: If possible, consolidate the following blocks
        if high_quality: 
            get_ss.command('\n'.join(
                    f'''
                    paste <( tabix -h {bgz_fname} {chrom} | \\
                            awk '{{print $1=$1":"$2":"$3":"$4, $2=${pval_col_idx}{awk_arg1}}}' | \\
                            sed -e 's/chr:pos:ref:alt/SNP/g' {sed_arg} ) \\
                          <( tabix -h {variant_manifest_bgz} {chrom} | \\
                            awk '{{ print $9 }}' ) | \\
                    awk '{{if({awk_arg2}) print $1,$2}}' > {get_ss[f"ofile_{chrom}"]}
                    '''
                    for chrom in chromosomes
                    )
            )
        else:
            get_ss.command('\n'.join(
                    f'''
                    tabix -h {bgz_fname} {chrom} | \\
                    awk '{{print $1=$1":"$2":"$3":"$4, $2=${pval_col_idx}{awk_arg1}}}' | \\
                    sed -e 's/chr:pos:ref:alt/SNP/g' {sed_arg} | \\
                    awk '{{if({awk_arg2}) print $1,$2}}' > {get_ss[f"ofile_{chrom}"]}
                    '''
                    for chrom in chromosomes
                    )
            )
                              
    ss_dict = {
            chrom: get_ss[f'ofile_{chrom}']
            for chrom in chromosomes
    }
        
    return ss_dict
    
def get_adj_betas(p, pop, not_pop, max_pops, pheno_key_dict, pheno_id, high_quality, 
                  hail_script):
    r'''
    Wrapper method for both PLINK clumping and SBayesR
    '''
    output_dir = (f'{ldprune_dir}/results{"_high_quality" if high_quality else ""}/'+
                  ('max_pops' if max_pops else f'{"not_" if not_pop else ""}{pop}')+
                  f'/{pheno_id}')

    clump_output_txt = f'{output_dir}/clump_results.txt' # PLINK clump output txt file
    clump_output_ht = f'{output_dir}/clump_results.ht' # PLINK clump output hail table
    
#    sbayesr_output_txt = f'{output_dir}/sbayesr_results-test.txt' # SBayesR output txt file
#    sbayesr_output_ht = f'{output_dir}/sbayesr_results-test.ht' # SBayesR output hail table
    
    overwrite = False
    
    clump_file_exists = hl.hadoop_is_file(f'{clump_output_ht}/_SUCCESS')
    if not clump_file_exists or overwrite:
        if clump_file_exists and overwrite:
            print('\n\nWARNING: Existing results will be overwritten for {pheno_id} in {output_dir}!\n')

        ss_dict = get_sumstats(p=p,
                               pop=pop, 
                               not_pop=not_pop,
                               max_pops=max_pops,
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
                   method='clump')
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
               hail_script, output_txt, output_ht, ss_dict, method):
    r'''
    Runs either PLINK clump (method = 'clump') or SBayesR (method = 'sbayesr')
    '''
    assert method in {'clump','sbayesr'}
    
    task_suffix = (f'{"not_" if not_pop else ""}{pop}' if not max_pops else 'max_pops')+f'-{pheno_id}'
    # TODO: if method = 'sbayesr' check if LD matrix has already been calculated
    
    tasks = []
    
    ref_subset = '-'.join(pheno_key_dict['pops'] if max_pops else ([p for p in POPS if p is not pop] if not_pop else [pop]))
    print(f'Using LD reference panel of {ref_subset}')
        
    ## run plink clumping
    for chrom, ss_chrom in ss_dict.items():
        ## read ref ld plink files 
        bfile = read_plink_input_group_chrom(p=p, 
                                             method=method,
                                             subset=ref_subset,
                                             chrom=chrom)
        
        get_betas = p.new_job(name=f'{method}_{task_suffix}_chr{chrom}')
        
        # TODO: change image to include GCTB if running SBayesR?
        get_betas.cpu(1) # plink clump cannot multithread
        
        get_betas.command('set -ex')
        
        if method == 'clump':
            get_betas.storage('5G') # default: 5G
#            clump_memory = -15*(chrom-1)+400 # Memory requested for PLINK clumping in MB. equation: -15*(chrom-1) + 500 is based on 400 MB for chr 1, 80 MB for chr 22
            clump_memory = 3.75 # in GB
            get_betas.memory(clump_memory) # default: 30G
            get_betas.command(f'head {ss_chrom}')
            get_betas.command(' '.join(['plink',
                                        '--bfile', str(bfile),
                                        '--memory',str(clump_memory*1000), # memory in MB
                                        '--threads','1', # explicitly set threads to 1
                                        '--clump', ss_chrom,
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

def main():    
    backend = hb.ServiceBackend(billing_project='ukb_diverse_pops',
                                   bucket='ukbb-diverse-temp-30day/nb-batch-tmp')
#    backend = batch.LocalBackend(tmp_dir='/tmp/batch/')

    high_quality = True
    not_pop = False
    max_pops = True
    
    paridx = 0
    parsplit = 1
    
    p = hb.batch.Batch(name=(f'clump{"-hq" if high_quality else ""}{"-max_pops" if max_pops else ""}'+
                             f'-{paridx}-{parsplit}'), 
                       backend=backend, default_image='gcr.io/ukbb-diversepops-neale/nbaya_plink:0.1',
                       default_storage='500Mi', default_cpu=8)
    
    ## download hail script to VM
    hail_script = p.read_input(f'{ldprune_dir}/scripts/python/plink_clump_hail.py')
    
    ## get phenotype manifest
    pheno_manifest = hl.import_table(f'{ldprune_dir}/phenotype_manifest.tsv.bgz',
                                     impute=True)
    pheno_manifest = pheno_manifest.to_pandas()
    
    if max_pops:
        pheno_list = get_pheno_list(pheno_manifest=pheno_manifest,
                                    pop=None,
                                    not_pop=False,
                                    max_pops=True)
    else:
        pheno_list = []
        for not_pop in [True, False]:
            for pop in POPS:
    
                pheno_list.append(get_pheno_list(pheno_manifest=pheno_manifest,
                                                 pop=pop,
                                                 not_pop=not_pop))
    
    idx_to_run = range(paridx, len(pheno_list), parsplit)
    pheno_list = [pheno_info for idx, pheno_info in enumerate(pheno_list) if idx in idx_to_run]
    print(f'Number of phens: {len(pheno_list)}')
    
    for trait_type, phenocode, pheno_sex, coding, modifier, pops, pop, not_pop, max_pops in pheno_list:
        pheno_key_dict = get_pheno_key_dict(trait_type, phenocode, pheno_sex, coding, modifier, pops)
        pheno_id = get_pheno_id(trait_type, phenocode, pheno_sex, coding, modifier)

        get_adj_betas(p=p,
                      pop=pop, 
                      not_pop=not_pop,
                      max_pops=max_pops,
                      pheno_key_dict=pheno_key_dict,
                      high_quality=high_quality,
                      pheno_id=pheno_id, 
                      hail_script=hail_script)
            
    p.run(open=True)
    
#    if type(backend)==batch.ServiceBackend:
#        print('running')
#    else:
#        p.run(verbose=True,
#              delete_scratch_on_exit=True)
    
    backend.close()


if __name__=="__main__":
    
    hl.init(default_reference='GRCh38', 
            spark_conf={'spark.hadoop.fs.gs.requester.pays.mode': 'AUTO', 
                        'spark.hadoop.fs.gs.requester.pays.project.id': 'ukbb-diversepops-neale'})
    main()
