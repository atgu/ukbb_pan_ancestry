#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 09:17:42 2020

Get PLINK subsets for clumping

@author: nbaya
"""

import hail as hl
import argparse
import hailtop.batch as hb
from ukbb_pan_ancestry.resources.genotypes import get_filtered_mt
from ukbb_pan_ancestry.resources.results import get_pheno_manifest_path
from ukbb_pan_ancestry import POPS, bucket

# MIN_CASES = 50
# MIN_CASES_ALL = 100
# MIN_CASES_EUR = 100

pop_dict = {
    'AFR': 6637, # dict with sample counts for each population 
    'AMR': 982,
    'CSA': 8876,
    'EAS': 2709,
    'EUR': 420542,
    'MID': 1599
            }

chroms = list(range(1,23))+['X']


def get_pops_list(args):
    r'''
    Generates list of population combinations. If `pops`=None, this will read 
    the phenotype manifest to get all population combinations.
    '''
    if args.pops is None:
        pheno_manifest = hl.import_table(get_pheno_manifest_path())
        pops_list_all = pheno_manifest.pops.collect()
        pops_list_all = sorted(list(set(pops_list_all)))
        pops_list_all = [p.split(',') for p in pops_list_all] # list of lists of strings
        idx = range(args.paridx, len(pops_list_all), args.parsplit)
        pops_list = [pops for i, pops in enumerate(pops_list_all) if i in idx]
    else:
        pops = sorted(args.pops.upper().split('-'))
        assert set(pops).issubset(POPS), f'Invalid populations: {set(pops).difference(POPS)}'
        pops_list = [pops] # list of list of strings
    print(f'''\n\npops: {'-'.join(pops) if len(pops_list)==1 else f"{len(pops_list)} of {len(pops_list_all)} combinations"}\n''')
    return pops_list


def get_bfile_chr_path(bfile_prefix, chrom):
    return f'{bfile_prefix}.chr{chrom}'


def get_mt_filtered_by_pops(pops: list,
                            chrom: str = 'all',
                            imputed: bool = True,
                            min_mac: int = 20,
                            entry_fields=('GP',),
                            filter_mac_instead_of_ac: bool = False):
    r'''
    Wraps `get_filtered_mt()` from ukbb_pan_ancestry.resources.genotypes
    This filters to samples from populations listed in `pops`.

    NOTE: If chrom='all', this loads all autosomes AND chrX.
    
    '''
    assert len(pops)>0 and set(pops).issubset(POPS)
    
    kwargs = {'pop': 'all' if len(pops)>1 else pops[0],
              'imputed': imputed,
              'min_mac': min_mac,
              'entry_fields': entry_fields,
              'filter_mac_instead_of_ac': filter_mac_instead_of_ac
              }
    
    mt = get_filtered_mt(chrom=chrom,  **kwargs) # in this case chrom='all' gets autosomes
    
    if chrom=='all':
        mt_x = get_filtered_mt(chrom='X', **kwargs)
        mt = mt.union_rows(mt_x)
    
    if len(pops)>1:
        mt = mt.filter_cols(hl.set(pops).contains(mt.pop))
        
    return mt


def get_pop_prop_dict(pop_dict: dict, pops: list) -> (dict, int):
    r'''
    Get population proportions in `pop_dict` for a list of populations `pops`
    '''
    tmp_pop_dict = {pop:n_pop for pop,n_pop in pop_dict.items() if pop in pops}
    n_total = sum(tmp_pop_dict.values())
    pop_prop_dict = {k: v/n_total for k,v in tmp_pop_dict.items()}
    return pop_prop_dict, n_total

def get_subset(mt_pop, pop_dict: dict, pops: list, n_max: int):
    r'''
    Get Hail table sample of max size = `n_max` for list of populations `pops`.
    '''
    pop_prop_dict, n_total = get_pop_prop_dict(pop_dict=pop_dict, 
                                               pops=pops)

    limiting_pop = min(pop_prop_dict, key=pop_prop_dict.get)
    n_sample = int(min(pop_dict[limiting_pop]/pop_prop_dict[limiting_pop], n_max))
    if n_sample != n_max:
        print(f'Using sample size of {n_sample} instead of {n_max} due to limiting population size in {limiting_pop}')
    print({k:v*n_sample for k,v in pop_prop_dict.items()})
        
    cols = mt_pop.cols()
    if len(pops)==1 and n_sample == pop_dict[pops[0]]: # if sampling a single population `pop` and n_sample is the same as the population's size. 
        ht_sample = cols
    else:
        cols = cols.annotate(tmp_rand = hl.rand_norm())
        cols = cols.order_by('tmp_rand')
        cols = cols.add_index(name = 'rand_idx')
        ht_sample = cols.filter(cols.rand_idx<n_sample)
        ht_sample = ht_sample.drop('tmp_rand','rand_idx')
    ht_sample = ht_sample.key_by('s')
    ht_sample = ht_sample.select('pop') # keyed by 's', thus the two remaining fields are 'pop' and 's'
    
    return ht_sample
    
def to_plink(pops: list, 
             subsets_dir, 
             mt, 
             ht_sample, 
             bfile_path,
             export_varid: bool = True,
             overwrite=False):
    r'''
    Exports matrix table to PLINK2 files
    NOTE: These files will need to split up by chromosome before plink_clump.py
    can be run. 
    '''
    assert 'GT' in mt.entry, "mt must have 'GT' as an entry field"
    assert mt.GT.dtype==hl.tcall, "entry field 'GT' must be of type `Call`"
    
    if not overwrite and all([hl.hadoop_exists(f'{bfile_path}.{suffix}') for suffix in ['bed','bim']]):
        print(f'\nPLINK .bed and .bim files already exist for {bfile_path}')
        print(bfile_path)
    else:
        print(f'Saving to bfile prefix {bfile_path}')
        mt_sample = mt.annotate_rows(varid = hl.str(mt.locus)+':'+mt.alleles[0]+':'+mt.alleles[1])
        mt_sample = mt_sample.filter_cols(hl.is_defined(ht_sample[mt_sample.s]))
        hl.export_plink(dataset = mt_sample, 
                        output = bfile_path, 
                        ind_id = mt_sample.s,
                        varid = mt_sample.varid) # varid used to be rsid

def export_varid(args):
    r'''
    Only used to check varids
    '''
    n_max = 5000
    subsets_dir = f'{bucket}/ld_prune/subsets_{round(n_max/1e3)}k' 
    
    mt = get_mt_filtered_by_pops(chrom='all', 
                                 pop='all',
                                 entry_fields=('GT',)) # default entry_fields will be 'GP', we need 'GT' for exporting to PLINK
    
    mt_sample = mt.annotate_rows(chrom = mt.locus.contig,
                                 pos = mt.locus.position,
                                 varid = hl.str(mt.locus)+':'+mt.alleles[0]+':'+mt.alleles[1])
    
    rows = mt_sample.rows()
    rows = rows.key_by()
    rows = rows.select('chrom','pos','varid')
    rows.export(f'{subsets_dir}/varid.txt',delimiter=' ')
    
def batch_split_by_chrom(args):
    r'''
    Splits bfiles by chromosome, for later use by plink_clump.py
    About $0.06 per population set
    '''
    
    hl.init(default_reference='GRCh38', 
            spark_conf={'spark.hadoop.fs.gs.requester.pays.mode': 'AUTO', 
                        'spark.hadoop.fs.gs.requester.pays.project.id': 'ukbb-diversepops-neale'})
    
    pops_list = get_pops_list(args)
    
    n_max = 5000 # maximum number of samples in subset (equal to final sample size if there are sufficient samples for each population)
    subsets_dir = f'{bucket}/ld_prune/subsets_{round(n_max/1e3)}k'
    
    backend = hb.ServiceBackend(billing_project='ukb_diverse_pops',
                                bucket='ukbb-diverse-temp-30day/nb-batch-tmp')
#    backend = batch.LocalBackend(tmp_dir='/tmp/batch/')
    
    b = hb.batch.Batch(name='split_by_chrom', backend=backend,
                       default_image='gcr.io/ukbb-diversepops-neale/nbaya_plink:0.1',
                       default_storage='30G', default_cpu=8)
        
    for pops in pops_list:
        pops_str = '-'.join(pops)
        bfile_prefix = f'{subsets_dir}/{pops_str}/{pops_str}'
        master_bfile_paths = [f'{bfile_prefix}.{suffix}' for suffix in ['bed','bim','fam']]
        master_fam_path = f'{bfile_prefix}.fam'
        bfile_chr_paths = [f'{get_bfile_chr_path(bfile_prefix, chrom)}.{suffix}' for chrom in chroms for suffix in ['bed','bim']]
        if not args.overwrite_plink and all(list(map(hl.hadoop_is_file, 
                                                      [master_fam_path]+bfile_chr_paths))):
            print(f'\nAll per-chrom PLINK files created for {pops_str}')
        else:
            if not all(map(hl.hadoop_is_file, master_bfile_paths)):
                print(f'\nWARNING: Insufficient files for {pops_str} to split into per-chrom bed/bim files, skipping\n')
                continue
            else:
                print(f'\n... Running bfile per-chrom split for {pops_str} ...')
                prefix = f'{subsets_dir}/{pops_str}/{pops_str}'
                bfile = b.read_input_group(
                    **{suffix:f'{prefix}.{suffix}' for suffix in ['bed','bim','fam']}
                    )
                split = b.new_job(name=f'split_by_chrom_{pops_str}')
                for chrom in chroms:
                    split.declare_resource_group(**{f'ofile_{chrom}':{'bed': '{root}.bed',
                                                                      'bim': '{root}.bim'}}) # exclude fam file to avoid redundancy
                    split.command(
                        f'''
                        plink \\
                        --bfile {bfile} \\
                        --chr {chrom} \\
                        --output-chr M \\
                        --make-bed \\
                        --out {split[f"ofile_{chrom}"]}
                        '''
                        )
                    # print(f"saving to {get_bfile_chr_path(bfile_prefix, chrom)}")
                    b.write_output(split[f'ofile_{chrom}'], get_bfile_chr_path(bfile_prefix, chrom))

    b.run(open=True)
    backend.close()

        
def main(args):
    
    hl.init(log='/tmp/hail.log')
    
    n_max = 5000 # maximum number of samples in subset (equal to final sample size if there are sufficient samples for each population)
    subsets_dir = f'{bucket}/ld_prune/subsets_{round(n_max/1e3)}k' 
    
    pops_list = get_pops_list(args)
    print(f'overwrite_plink: {args.overwrite_plink}')

    for pops in pops_list:
        pops_str = '-'.join(pops)
        ht_sample_path = f'{subsets_dir}/{pops_str}/{pops_str}.ht'
        bfile_prefix = f'{subsets_dir}/{pops_str}/{pops_str}'
        
        master_bfile_paths = [f'{bfile_prefix}.{suffix}' for suffix in ['bed','bim','fam']]
        
        if not args.overwrite_plink and all(map(hl.hadoop_is_file, 
                                                [f'{ht_sample_path}/_SUCCESS']+master_bfile_paths)):
            continue
        else:
            print(f'\n... Starting PLINK exports for {pops_str} ...')
            mt_pop = get_mt_filtered_by_pops(pops=pops,
                                             chrom='all',  # chrom='all' includes autosomes and chrX
                                             entry_fields=('GT',)) # default entry_fields will be 'GP', we need 'GT' for exporting to PLINK
            if hl.hadoop_is_file(f'{ht_sample_path}/_SUCCESS'):
                ht_sample = hl.read_table(ht_sample_path)
                ht_sample_ct = ht_sample.count()
                print(f'... Subset ht already exists for pops={pops_str} ...')
                print(f'\nSubset ht sample ct: {ht_sample_ct}\n\n')
            else:
            
                print(f'\n\n... Getting sample subset ({pops_str}) ...\n')
                
                ht_sample = get_subset(mt_pop = mt_pop,
                                        pop_dict = pop_dict, 
                                        pops = pops, 
                                        n_max = n_max)
                
                ht_sample_ct = ht_sample.count()
                print(f'\n\nht_sample_ct: {ht_sample_ct}\n\n')
                ht_sample = ht_sample.checkpoint(ht_sample_path)
            
            print(f'... Exporting to PLINK ({pops_str}) ...')
            to_plink(pops = pops,
                      subsets_dir=subsets_dir,
                      mt = mt_pop,
                      ht_sample = ht_sample,
                      bfile_path = bfile_prefix,
                      overwrite=args.overwrite_plink)
        

if __name__=='__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--pops', default=None, type=str, help='population to use')
    parser.add_argument('--overwrite_plink', default=False, action='store_true', help='whether to overwrite existing PLINK files')
    parser.add_argument('--export_varid', default=False, action='store_true', help='export varids')
    parser.add_argument('--batch_split_by_chrom', default=False, action='store_true', help='Whether to split PLINK files into per-chrom files')
    parser.add_argument('--parsplit', type=int, default=1, help="number of parallel batches to split pop combinations into")
    parser.add_argument('--paridx', type=int, default=0, help="which of the parallel batches to run (zero-indexed)")
    args = parser.parse_args()
    
    if args.export_varid:
        export_varid(args=args)
    if args.batch_split_by_chrom:
        batch_split_by_chrom(args)
    else:
        main(args=args)