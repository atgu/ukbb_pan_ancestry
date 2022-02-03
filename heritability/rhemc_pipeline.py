__author__ = 'Rahul Gupta'

import hail as hl

# hl.init(spark_conf={'spark.hadoop.fs.gs.requester.pays.mode': 'AUTO',
#                     'spark.hadoop.fs.gs.requester.pays.project.id': 'ukbb-diversepops-neale'})

import hailtop.batch as hb
import numpy as np
import math
import logging, time
import sys
import re
import argparse

from os.path import basename
from typing import Union
from google.cloud import storage

from ukbb_pan_ancestry.get_plink_subsets import get_mt_filtered_by_pops
from ukbb_pan_ancestry.resources.genotypes import get_ukb_af_ht_path
from ukbb_pan_ancestry.resources.ld import get_ld_score_ht_path
from ukbb_pan_ancestry.resources.phenotypes import get_ukb_pheno_mt_path
from ukb_common.resources.generic import PHENO_KEY_FIELDS
from ukbb_pan_ancestry.resources.results import get_variant_results_path

# paths
bucket = 'ukb-diverse-pops'
fusebucket = 'rgupta-pcgc-mount'
loc = 'rg-pcgc'
path_geno = f'gs://{bucket}/{loc}/genos/'
path_pheno = f'gs://{bucket}/{loc}/phenos/'
path_covar = f'gs://{bucket}/{loc}/covar/'
path_annot = f'gs://{bucket}/{loc}/annot/'
path_results = f'gs://{bucket}/{loc}/results/'

pheno_script_adr = f'gs://{bucket}/{loc}/code/rhemc_load_pheno.py'
cat_script_adr = f'gs://{bucket}/{loc}/code/concat_tables.py'
parsing_script_adr = f'gs://{bucket}/{loc}/code/rhe_mc_parse_results.py'

SAIGE_PHENO_LOC = f'gs://{bucket}/results/pheno_export_data/'

MT_TEMP_BUCKET = f'gs://ukbb-diverse-temp-30day/{loc}/mt/'
HT_TEMP_BUCKET = f'gs://ukbb-diverse-temp-30day/{loc}/ht/'

# other environment variables
IMAGE = 'gcr.io/ukbb-diversepops-neale/rgupta_rhe_mc'
MHC = hl.parse_locus_interval('6:25M-35M') # MHC locus
PHWE = 1e-7 # p HWE cutoff
RANDOMPREF = lambda idx: f'random-stdnormal-{str(idx)}-both_sexes'

# MAF ranges for bins. Recall that LD scores bins are created via quantiles.
MAFRANGE_5 = [(0.0, 0.1),(0.1,0.2),(0.2,0.3),(0.3,0.4),(0.4,1.0)]
MAFRANGE_2 = [(0.0,0.05),(0.05,1.0)]


def check_indiv_files(ancestries, args, nbins):
    """ Check if individual level files exist. This means genotype, covariates, 
    and annotations. If files are not found, creates them from source data.
    This *does not* verify if phenotype files are present.

    Parameters
    ----------
    ancestries: `list`
    List of ancestries to include.
    """
    logging.info('Checking for genotype and annotation files...')
    # genotype and annotations
    geno_names = get_geno_split_names(ancestries, dictout=False)
    tf_geno = [hl.hadoop_exists(path_geno + x) for x in geno_names]
    annot_names = get_annot_split_names(ancestries, dictout=False, n_annot=nbins)
    tf_annot = [hl.hadoop_exists(path_annot + x) for x in annot_names]
    if not all(tf_geno) or not all(tf_annot):
        if not all(tf_geno):
            genomiss = ', '.join([nm for nm, tf in zip(geno_names, tf_geno) if not tf])
            logging.info('The following genotype files were not found: ' + genomiss)
        if not all(tf_annot):
            annotmiss = ', '.join([nm for nm, tf in zip(annot_names, tf_annot) if not tf])
            logging.info('The following annotation files were not found: ' + annotmiss)
        # if args.verbose:
        #     print('Full set of genotype and/or annotation files not found. Generating...')
        generate_geno_annot_split(path_geno, path_annot, ancestries, args, nbins)
        logging.info('Genotype and annotation files generated.')
        # if args.verbose:
        #     print('Genotype and annotation files generated.')
    else:
        logging.info('All genotype and annotation files found.')
        # if args.verbose:
        #     print('Genotype and annotation files found.')

    # Ensure that .fam files are converted into Hail Tables and that the phenotype table exists
    logging.info('Validating/setting up Hail tables from .fam files for future use...')
    _ = get_famfiles(ancestries, checkpoint=True, override_check=False)
    logging.info('Validating/setting up phenotype MatrixTable for future use...')
    _ = _read_pheno_data(checkpoint=True)

    # covariate
    logging.info('Checking for covariate files...')
    covar_names = get_covar_split_names(ancestries, dictout=False, sex_specific=args.sex_specific)
    tf_covar = [hl.hadoop_exists(path_covar + x) for x in covar_names]
    if not all(tf_covar):
        covarmiss = ', '.join([nm for nm, tf in zip(covar_names, tf_covar) if not tf])
        logging.info('The following covariate files were not found: ' + covarmiss)
        generate_covar_split(path_covar, ancestries, args.checkpoint, args.sex_specific)

    else:
        logging.info('All covariate files found.')
    
    logging.info('All covariate, annotation, and genotype files generated.')


def gscopy(source, dest):
    """Copies a blob from one bucket to another with a new name."""
    match_source = re.search(r'gs://([A-Za-z1-9\-]{1,})/(.+)',source)
    match_dest = re.search(r'gs://([A-Za-z1-9\-]{1,})/(.+)',dest)
    src_bucket = match_source.group(1)
    src_blob = match_source.group(2)
    dest_bucket = match_dest.group(1)
    dest_blob = match_dest.group(2)

    storage_client = storage.Client()

    source_bucket = storage_client.bucket(src_bucket)
    source_blob = source_bucket.blob(src_blob)
    destination_bucket = storage_client.bucket(dest_bucket)

    _ = source_bucket.copy_blob(source_blob, destination_bucket, dest_blob)


def prepare_fuse_bucket(ancestries, nbins, sex_specific):
    """ Copies genotype, covariate, and annotation files over to the fuse bucket.
    """
    logging.info(f'Preparing bucket {fusebucket} by adding covariate, genotype, and annotation files...')
    covar_files = get_covar_split_names(ancestries, True, sex_specific)
    geno_files = get_geno_split_names(ancestries, True, False)
    annot_files = get_annot_split_names(ancestries, True, n_annot=nbins)
    log_copy = lambda src, dest: logging.info(f'Copying {src} to {dest}...')

    for anc in ancestries:
        # covariate
        source_covar = path_covar + covar_files[anc]
        dest_covar = f'gs://{fusebucket}/covar/{covar_files[anc]}'
        if not hl.hadoop_exists(dest_covar):
            log_copy(source_covar, dest_covar)
            gscopy(source_covar, dest_covar)

        # genotype
        source_geno = [path_geno + x for x in geno_files[anc]]
        dest_geno = [f'gs://{fusebucket}/genos/{x}' for x in geno_files[anc]]
        for idx in range(0, len(dest_geno)):
            if not hl.hadoop_exists(dest_geno[idx]):
                log_copy(source_geno, dest_geno)
                gscopy(source_geno[idx], dest_geno[idx])

        # annot
        source_annot = path_annot + annot_files[anc]
        dest_annot = f'gs://{fusebucket}/annot/{annot_files[anc]}'
        if not hl.hadoop_exists(dest_annot):
            log_copy(source_annot, dest_annot)
            gscopy(source_annot, dest_annot)

    logging.info(f'All files successfully transferred or verified.')


def make_quantiles(expr, nq=5, approx=False, ht=None):
    """ Makes quantiles from `Expression` in expr.

    Parameters
    ----------
    expr: `Expression`
    To be broken up into quantiles.

    nq: `int32`
    Number of quantiles to create.

    approx: `bool`
    If appox, will use the Hail `approx_quantiles` function. If not, will
    localize the expression to identify quantiles locally. If True, ht must be specified.
    
    Returns
    -------
    `list` of `tuples`
    """
    pr_vec = [i/nq for i in range(0,nq+1)]
    if approx:
        if ht is None: raise TypeError('ht must be specified if using the Hail quantile approximation.')
        lower, upper, quants = ht.aggregate((hl.agg.min(expr), 
                                             hl.agg.max(expr),
                                             hl.agg.approx_quantiles(expr, pr_vec, 1000)))
    else:
        value_vector = [i for i in expr.collect() if i is not None]
        lower = np.min(value_vector)
        upper = np.max(value_vector)
        quants = list(np.quantile(value_vector, pr_vec))
    quants[len(quants)-1] = upper + 1 # expand outwards to ensure inclusivity later
    quants[0] = lower - 1 # expand outwards to ensure inclusivity later
    return list(zip(quants[0:nq], quants[1:(nq+1)]))


def discretize(expr, list_of_tuples):
    """ Takes an expression and cuts up results into bins specified by list_of_tuples.
    Anything out of range is given an assignment of missing. 1-indexed.

    NOTE This does not check that intervals are non-overlapping. Must ensure this or else
    this function will return the first instance of this.

    Parameters
    ----------
    expr: `Expression`
    Hail expression to convert into discrete bins.

    list_of_tuples: `list` of `tuples`
    Contains boundary definitions of bins. Right side is assumed inclusive, left side exclusive.

    Returns
    -------
    `IntegerExpression`
    """
    def within(ele, tupe):
        lo, hi = tupe
        return ((ele > lo) & (ele <= hi))
    
    converted = hl.find(lambda y: within(expr, y[1]), list(enumerate(list_of_tuples)))
    return converted[0]+1


def explode_bins(ht, cols_bins, verbose):
    """ Explodes columns representing bin assignments into wide format. The columns
    in cols_bins will be removed and new columns will be labeled accordingly:

    >>> cols_bins = [maf_bin, LD_bin]
    Produces new columns: `maf_bin_1_LD_bin_1; maf_bin_2_LD_bin_1; ....; maf_bin_5_LD_bin_5`

    New columns will be of type int32, containing 1 if that record belongs to that bin
    and otherwise containing 0. If the bin assignment is missing, all columns will list 0.

    NOTE Requires that all fields in `cols_bins` are `int32` and are 1-indexed consecutive
    integers. None or missing values are acceptable.

    Parameters
    ----------
    ht: `HailTable`

    cols_bins: `list` of `str`
    Contains field names found in ht which are to be exploded.

    verbose: `bool`
    If True, will output information on per-bin counts

    Returns
    -------
    `HailTable`
    """
    # check types
    if not all([ht[x].dtype == hl.dtype('int32') for x in cols_bins]):
        raise TypeError('All fields in cols_bins must have type int32.')
    # ensure sequential 1-indexed integers
    vals_map = {}
    for ele in cols_bins:
        n_count = ht.aggregate(hl.agg.counter(ht[ele]))
        if verbose:
            print('Bin distribution across variants for ' + ele + ':')
            print(n_count)
        vals_vec = [i for i in list(n_count.keys()) if i is not None]
        vals_vec.sort()
        tf = vals_vec == [i for i in range(1,len(vals_vec)+1)]
        if not tf:
            raise ValueError('Field ' + ele + ' is not consecutive 1-indexed. ' + \
                             'Values found were: ' + str(vals_vec))
        vals_map.update({ele: vals_vec})
    
    # transform to wide format
    holder = {}
    for ele in cols_bins:
        if len(holder) == 0:
            holder = {ele + '_' + str(x): [x] for x in vals_map[ele]}
        else:
            holder = {item + '_' + ele + '_' + str(x): arr + [x] 
                      for item, arr in holder.items() for x in vals_map[ele]}
    
    # create new columns
    col_holder = {}
    for k,v in holder.items():
        conditions = [ht[ele] == val for ele, val in zip(cols_bins,holder[k])]
        final_cond = conditions[0]
        for i in range(1, len(conditions)):
            final_cond = final_cond & conditions[i]
        col_holder.update({k: hl.if_else(final_cond, 1, 0)})
    ht = ht.annotate(**col_holder)
    
    # remove original columns
    ht = ht.drop(*cols_bins)
    
    # return
    return ht


def generate_geno_annot_split(path_geno, path_annot, ancestries, args, nbins):
    """ Generates genotype and annotation files. These will be split by ancestry.
    This function guarentees that each ancestry will have the same set of SNPs, and that
    the annotations (which are different for each ancestry) will contain the same SNPs
    in the same order as those listed in the PLINK data.

    Parameters
    ----------
    path_geno: `str`
    Path to genotype files.

    path_annot: `str`
    Path to annotation files.

    ancestries: `list`
    List of ancestries to include

    args: arguments to pipeline
    """
    logging.info('Generating genotype files...')
    # obtain genotype data
    geno_full_path = MT_TEMP_BUCKET + 'initial_genotypes_all_ancestry.mt'
    if hl.hadoop_exists(geno_full_path):
        logging.info('Genotype (prefiltered by pops and withdrawal) MatrixTable found and loaded.')
        mt = hl.read_matrix_table(geno_full_path)
    else:
        logging.info('Creating genotype (prefiltered by pops and withdrawal) MatrixTable...')
        mt = get_mt_filtered_by_pops(pops=ancestries, entry_fields=('GT',))
        if args.checkpoint:
            mt = mt.checkpoint(geno_full_path, overwrite=True)

    if args.verbose:
        print('Imported per-population genotype Ns:')
        print_pop_Ns(mt)

    MAFRANGE = MAFRANGE_2 if args.maf_bins_2 else MAFRANGE_5
        
    geno_filtered_path = MT_TEMP_BUCKET + 'filtered_genotypes_all_ancestry.mt'
    if hl.hadoop_exists(geno_filtered_path):
        logging.info('Genotype MatrixTable, filtered, loaded.')
        mt_filt = hl.read_matrix_table(geno_filtered_path)
    else:
        logging.info('Filtering genotype MatrixTable...')
        # import variant level data
        af_ht = hl.read_table(get_ukb_af_ht_path())

        # filter MAF > cutoff (in all populations) and is defined (all populations)
        af_ht_f = af_ht.filter(hl.all(lambda x: hl.is_defined(af_ht.af[x]), 
                                      hl.literal(ancestries)))
        af_ht_f = af_ht_f.filter(hl.all(lambda x: (af_ht.af[x] >= args.maf) & \
                                                  (af_ht.af[x] <= (1-args.maf)), 
                                      hl.literal(ancestries)))
        mt_maf = mt.filter_rows(hl.is_defined(af_ht_f[mt.row_key]))

        # remove relateds
        mt_nonrel = mt_maf.filter_cols(~mt_maf.related)

        # compute phwe, remove those with p < 1e-7
        mt_nonrel_hwe = mt_nonrel.annotate_rows(**{'hwe_' + anc.lower(): 
                                                    hl.agg.filter(mt_nonrel.pop == anc, 
                                                                  hl.agg.hardy_weinberg_test(mt_nonrel.GT)) 
                                                    for anc in ancestries})
        ancestries_tf = [mt_nonrel_hwe['hwe_' + anc.lower()].p_value >= PHWE for anc in ancestries]
        mt_nonrel_hwe = mt_nonrel_hwe.filter_rows(hl.all(lambda x: x, ancestries_tf))

        # remove MHC region
        mt_filt = hl.filter_intervals(mt_nonrel_hwe, [MHC], keep=False)

        if args.checkpoint:
            mt_filt = mt_filt.checkpoint(geno_filtered_path, overwrite=True)

        logging.info('Filtering complete.')

    if args.verbose:
        print('Post-filtering per-population genotype Ns:')
        _ = print_pop_Ns(mt_filt)

    # output results as plink files
    output_files = get_geno_split_names(ancestries, dictout=True, plinkprefix=True)
    output_files_plink = get_geno_split_names(ancestries, dictout=True, plinkprefix=False)
    for anc in ancestries:
        plinkfiles = output_files_plink[anc]
        if all([hl.hadoop_exists(path_geno + x) for x in plinkfiles]):
            print('PLINK files for ' + anc + ' already exist. Skipping generation...')
        else:
            logging.info('Saving PLINK files for ' + anc + '...')
            mt_filt_anc = mt_filt.filter_cols(mt_filt.pop == anc)
            hl.export_plink(mt_filt_anc, path_geno+output_files[anc], ind_id = mt_filt_anc.s)

    # import LD score information & arrange to have the same order as genotype data
    # then produce MAF bins and LD score bins
    logging.info('Constructing annotation files...')
    snps_out = mt_filt.rows().select('varid')
    output_files_annot = get_annot_split_names(ancestries, dictout=True, n_annot=nbins)
    output_files_annot_noextn = get_annot_split_names(ancestries, dictout=True, n_annot=nbins, suffix_incl=False)
    for anc in ancestries:
        ht_anc = hl.read_table(get_ld_score_ht_path(pop=anc))
        ht_anc_expr = ht_anc[snps_out.row_key]
        this_tab = snps_out.annotate(ld_score = ht_anc_expr.ld_score,
                                     af = ht_anc_expr.AF)
        this_tab = this_tab.annotate(maf = 0.5 - hl.abs(0.5 - this_tab.af))
        LD_bin_anc = make_quantiles(this_tab.ld_score, nq=args.num_ld_bins, 
                                    approx=args.approx_quantiles, ht=this_tab)
        this_tab = this_tab.select(maf_bin = discretize(this_tab.maf, MAFRANGE),
                                   LD_bin = discretize(this_tab.ld_score, LD_bin_anc))
        exploded_tab = explode_bins(this_tab, ['maf_bin','LD_bin'], args.verbose)

        # back up annotations as a hail table
        exploded_tab.write(path_annot+'ht/'+output_files_annot_noextn[anc] + '.ht',
                           overwrite=True)

        # output
        _write_flat_without_key(exploded_tab, path_annot+output_files_annot[anc], 
                                delimiter=' ', header=False)
    logging.info('Annotation files created.')


def hl_combine_str(*expressions, sep="_") -> hl.StringExpression:
    """ This function combines StringExpressions. Notably, if a string is missing, 
    nothing is included.
    """
    strings = hl.array(list(expressions))
    return hl.delimit(strings.filter(lambda x: hl.is_defined(x) & (hl.len(x) > 0)), sep)

    # expressions_replaced = [hl.if_else(hl.is_missing(this_expr), "", this_expr) 
    #                             for this_expr in expressions]
    # final_expr = []
    # for idx,this_expr in enumerate(expressions_replaced):
    #     if idx == 0:
    #         final_expr = this_expr
    #     else:
    #         tf_1 = hl.len(final_expr) == 0
    #         tf_2 = hl.len(this_expr) == 0
    #         final_expr = (hl.case().when(tf_1 & tf_2, final_expr)
    #                                .when(tf_1, this_expr)
    #                                .when(tf_2, final_expr)
    #                                .default(final_expr + sep + this_expr))
    #         # final_expr = hl.if_else(tf_1 & tf_2, final_expr,
    #         #                         hl.if_else(tf_1, this_expr,
    #         #                                    hl.if_else(tf_2, final_expr, 
    #         #                                               final_expr + sep + this_expr)))
    # return final_expr


def _generate_map_anc_pheno(path_prefix, ancestries):
    dictout = {}
    for anc in ancestries:
        thisdir = path_prefix + anc + '/'
        if hl.hadoop_is_dir(thisdir):
            dictout.update({anc: [basename(x['path']) for x in hl.hadoop_ls(thisdir)]})
        else:
            dictout.update({anc: []})
    return dictout


def list_precomputed_pheno_files(ancestries):
    """ Obtain the set of phenotype files that have been pre-formatted for RHEmc.

    Parameters
    ----------
    ancestries: `list`

    Returns
    -------
    `dict`
    """
    return _generate_map_anc_pheno(path_pheno, ancestries)


def list_saige_pheno_files(ancestries):
    """ Obtain the set of phenotype files that were generated for the Saige run.

    Parameters
    ----------
    ancestries: `list`

    Returns
    -------
    `dict`
    """
    return _generate_map_anc_pheno(SAIGE_PHENO_LOC, ancestries)


def list_completed_rhemc_logs():
    """ Obtain the list of logs available from previous runs.
    """
    thisdir = path_results + 'raw/'
    if hl.hadoop_is_dir(thisdir):
        listout = [basename(x['path']) for x in hl.hadoop_ls(thisdir)]
    else:
        listout = []
    
    return listout



def construct_phenotype_id(tab: Union[hl.Table, hl.MatrixTable]):
    """ Returns a StringExpression representing a combination of phenotype keys.

    Parameters
    ----------
    tab: Table or MatrixTable
    Must contain all relevant phenotype keys for combination.

    Returns
    -------
    StringExpression
    """
    return hl_combine_str(*[tab[x] for x in PHENO_KEY_FIELDS], sep='-')
    #return hl.delimit([tab[x] for x in PHENO_KEY_FIELDS], '_')


def get_famfiles(ancestries, checkpoint=False, override_check=False):
    """Obtains fam files in HailTable format, trimmed to just include FID and IID.
    Throws error if not found. If fam files do not exist, run generate_geno_annot_split.

    Parameters
    ----------
    ancestries: `list`

    checkpoint: `bool`

    override_check: `bool`
    If True, will assume that plink files and the fam files already exist.

    Returns
    -------
    `dict` indexed by ancestries.
    """
    geno_files = get_geno_split_names(ancestries, dictout=True, plinkprefix=True)
    if not override_check:
        tf_geno = {anc: hl.hadoop_exists(path_geno + geno_files[anc] + '.fam') for anc in ancestries}
        if not all(tf_geno.values()):
            raise ValueError('All genotype plink files must exist. The following ancestries are missing: ' + \
                            ','.join([anc for anc in ancestries if not tf_geno[anc]]))
    
    # load fam files and make backbone for covariate and phenotype files
    famfiles = {}
    for anc in ancestries:
        thisloc = path_geno + 'ht/' + anc + '_famfile.ht'
        if override_check or hl.hadoop_exists(thisloc):
            famfiles.update({anc: hl.read_table(thisloc)})
        else:
            famtab = hl.import_table(path_geno + geno_files[anc] + '.fam', impute=True, no_header=True)
            famtab = famtab.select(FID=famtab.f0, IID=famtab.f1
                          ).add_index(name='idx'
                          ).key_by('idx'
                          ).repartition(100
                          ).order_by('idx'
                          ).drop('idx')
            if checkpoint:
                famtab = famtab.checkpoint(thisloc)
            famfiles.update({anc: famtab})
    
    return famfiles


def generate_covar_split(path_covar, ancestries, checkpoint=False, sex_specific=False):
    """This will construct and output covariate files. This function
    guarentees that the same individuals (in the same order) are included per person
    as in the genotype files (particularly the .fam file).
    NOTE: expects that genotype PLINK files have already been created. Will enforce indiviudal order
    based on this.

    Parameters
    ----------
    path_cvoar: `str`
    Path to covariate files.

    path_geno: `str`
    Path to genotype files.
    
    ancestries: `list`
    List of ancestries to include
    """
    if sex_specific:
        logging.info(f'Generating covariate files without sex covariates...')
        other_covars = ['age', 'age2']
    else:
        logging.info(f'Generating covariate files...')
        other_covars = ['age', 'sex', 'age_sex', 'age2', 'age2_sex']
    famfiles = get_famfiles(ancestries, checkpoint, override_check=False)
    
    # construct covariates and output
    covars = _read_pheno_rowtable()
    covariate_names = ['PC' + str(x) for x in range(1,21)] + other_covars
    covars = covars.select(*covariate_names)
    covar_dirs = get_covar_split_names(ancestries, dictout=True, sex_specific=sex_specific)
    for anc in ancestries:
        famfile_this = famfiles[anc]
        famfile_this = famfile_this.add_index('idx').key_by('IID')
        famfile_annot = famfile_this.join(covars, how='left')
        famfile_annot = famfile_annot.key_by('idx').order_by('idx').drop('idx')
        famfile_annot = famfile_annot.select(*(['FID','IID'] + covariate_names))
        famfile_annot.export(output=path_covar + covar_dirs[anc], header=True, delimiter='\t')


def get_geno_split_names(ancestries, dictout, plinkprefix=False):
    """ Returns names of genotype files.

    Parameters
    ----------
    ancestries: `list`
    List of ancestries to include

    dictout: `bool`
    If the output should be a dictionary indexed on ancestries.

    plinkprefix: `bool`
    If True, will return the plink file prefix. If False, will return a list of files (bim, bed, fam).
    """
    file_base = 'pananc_31063_pcgc_genos'
    plink_set = ['.bim', '.bed', '.fam']
    if dictout:
        if plinkprefix:
            out = {anc: anc + '_' + file_base for anc in ancestries}
        else:
            out = {anc: [anc + '_' + file_base + suff for suff in plink_set] for anc in ancestries}
    else:
        if plinkprefix:
            out = [anc + '_' + file_base for anc in ancestries]
        else:
            out = [anc + '_' + file_base + suff for anc in ancestries for suff in plink_set]
    
    return out


def get_pheno_split_names(ancestries, dictout, phenotype_id, enable_suffix=True):
    """ Returns names of phenotype files. The file name format is:
    ANCESTRY_phenotypeID_pananc_31063_pcgc_pheno.tsv.bgz

    Parameters
    ----------
    ancestries: `list`
    List of ancestries to include

    dictout: `bool`
    If the output should be a dictionary indexed on ancestries.

    phenotype_id: `str`
    The ID of the phenotype file. See construct_phenotype_id for how these are constructed.
    """
    file_base = 'pananc_31063_pcgc_pheno'
    suffix = '.tsv.bgz'
    pheno_suffix = phenotype_id + '_' + file_base + suffix if enable_suffix else phenotype_id + '_' + file_base
    if dictout:
        out = {anc: anc + '/' + pheno_suffix for anc in ancestries}
    else:
        out = [anc + '/' + pheno_suffix for anc in ancestries]

    return out


def get_pheno_filename(phenotype_id, enable_suffix=True):
    """ Returns names of phenotype files. These are the full file containing all individuals.
    These are not ancestry specific.

    Parameters
    ----------
    phenotype_id: `str`
    The ID of the phenotype file. See construct_phenotype_id for how these are constructed.
    """
    file_base = f'pananc_31063_pcgc_pheno'
    suffix = '.tsv.bgz'
    return phenotype_id + '_' + file_base + suffix if enable_suffix else phenotype_id + '_' + file_base


def get_covar_split_names(ancestries, dictout, sex_specific):
    """ Returns names of covariate files.

    Parameters
    ----------
    ancestries: `list`
    List of ancestries to include

    dictout: `bool`
    If the output should be a dictionary indexed on ancestries.

    sex_specific: `bool`
    If sex-specific covariates should be ELIMINATED (because the phenotype is sex-stratified).
    """
    file_base = 'pananc_31063_pcgc_covariates'
    file_base = file_base + ('_nosexcovar' if sex_specific else '')
    suffix = '.tsv'
    if dictout:
        out = {anc: anc + '_' + file_base + suffix for anc in ancestries}
    else:
        out = [anc + '_' + file_base + suffix for anc in ancestries]

    return out


def get_annot_split_names(ancestries, dictout, n_annot=25, suffix_incl=True):
    """ Returns names of annotation files.

    Parameters
    ----------
    ancestries: `list`
    List of ancestries to include

    dictout: `bool`
    If the output should be a dictionary indexed on ancestries.
    """
    file_base = f'pananc_31063_pcgc_{str(n_annot)}LDMS_annotations'
    suffix = '.txt.bgz' if suffix_incl else ''
    if dictout:
        out = {anc: anc + '_' + file_base + suffix for anc in ancestries}
    else:
        out = [anc + '_' + file_base + suffix for anc in ancestries]

    return out
    

def get_result_log_split_names(ancestries, phenotype_id, dictout, suffix=''):
    """ Returns names of the raw log to output from the RHEmc run.

    Parameters
    ----------
    ancestries: `list`
    List of ancestries to include

    dictout: `bool`
    If the output should be a dictionary indexed on ancestries.

    suffix: `str`
    Appended to the end of the filename before the extension.
    """
    file_base = 'result'
    suffix = f'{suffix}.log'
    if dictout:
        out = {anc: anc + '_' + phenotype_id + '_' + file_base + suffix for anc in ancestries}
    else:
        out = [anc + '_' + phenotype_id + '_' + file_base + suffix for anc in ancestries]

    return out


def parse_ancestries(args):
    return args.ancestries.split(',')


def construct_iter_suffix(iter):
    return '-iter' + str(iter) if iter is not None else ''


def print_pop_Ns(mt):
    """ Enumerates the number of individuals per ancestry represented in the genotype table.

    Parameters
    ----------
    mt: :obj: `MatrixTable`

    Returns
    -------
    None
    """
    dict_anc = mt.aggregate_cols(hl.agg.counter(mt.pop))
    _ = [print(k + ': ' + str(v)) for k,v in dict_anc.items()]


def convert_pheno_id_to_potential_saige(pheno_id: Union[hl.StringExpression, str]):
    """ Convert phenotype ID to string for matching to Saige.
    """
    sexes = ['both_sexes', 'females', 'males']
    exclude = 'phecode'
    split_pheno = pheno_id.split('-')
    if type(pheno_id) == str:
        if exclude in split_pheno:
            newstr = pheno_id
        else:
            newstr = '-'.join([x for x in split_pheno if x not in sexes])
    else:
        hljoin = lambda x: hl.str('-').join(x)
        newstr = hl.if_else(split_pheno.contains(hl.literal(exclude)), 
                            hljoin(split_pheno),
                            hljoin(hl.filter(lambda x: ~hl.literal(sexes).contains(x), split_pheno)))
    return newstr + '.tsv'


def run_phenotype_job(b, phenotype_id, ancestries, use_saige, checkpoint, 
                      phenoscript, n_threads, random):
    """ Runs phenotype file creation jobs.
    """
    j = b.new_job(name=phenotype_id+'_create_pheno')
    j.image(IMAGE)
    j.cpu(n_threads)
    filename = get_pheno_filename(phenotype_id, enable_suffix=False)
    filename_map = get_pheno_split_names(ancestries, dictout=True, phenotype_id=phenotype_id, enable_suffix=True)
    filename_compat = get_pheno_filename(compatiblify_phenotype_id(phenotype_id), enable_suffix=False)
    filename_map_compat = get_pheno_split_names(ancestries, dictout=True, 
                                                phenotype_id=compatiblify_phenotype_id(phenotype_id), 
                                                enable_suffix=True)
    ancestry = ','.join(ancestries)
    checkpoint_val = '--checkpoint' if checkpoint else ''
    saige_pull = '--pull-from-saige' if use_saige else ''
    random_val = '--random' if random else ''
    #write this file get_pheno_filename(args.phenotype_id, enable_suffix=False)
    anc_array =  '( ' + ', '.join(ancestries) + ' )'
    map_dependencies = {anc: j[anc] for anc in ancestries}
    command = f"""
        ancestry_array={anc_array}
        for i in '${{ancestry_array[@]}}'
        do
            mkdir $i
        done
        python3 {phenoscript} --phenotype-id '{phenotype_id}' --ancestries {ancestry} \
                              --logging --override-check {checkpoint_val} {saige_pull} {random_val}
        cp '{filename_compat + '.log'}' {j.log}
              """
    for anc in ancestries:
        command+=f"""
            cp '{filename_map_compat[anc]}' {map_dependencies[anc]}
        """
    j.command(command)
    b.write_output(j.log, path_pheno + 'log/' + filename + '.log')
    for anc in ancestries:
        b.write_output(map_dependencies[anc], path_pheno + filename_map[anc])
    
    return b, j, map_dependencies


def run_rhemc_job(b, phenotype_id, ancestry, phenotype_file, read_previous: bool, parser, 
                 use_fuse: bool, memsize: int, storage, jackknife_blocks: int, 
                 random_vectors: int, nbins: int, suffix: str, iter, sex_specific):
    # Localize all files (and pass phenotype files in from previous step)
    # This will require GCSFUSE for the annotation, covariate, and genotype files
    # Run RHEmc 
    # Collect results into a table
    iter_suffix = construct_iter_suffix(iter)
    j = b.new_job(name=ancestry + '_' + phenotype_id + '_RHEmc' + iter_suffix)
    j.image(IMAGE)
    log_out = get_result_log_split_names([ancestry], phenotype_id, dictout=True, suffix=suffix)[ancestry] + iter_suffix
    if read_previous:
        log_file = b.read_input(path_results + 'raw/' + log_out)       
        command = f"""
            python3 {parser} \
                --file {log_file} \
                --pheno '{phenotype_id}{iter_suffix}' \
                --ancestry {ancestry} \
                --out {j.parsed_out}
        """
        j.command(command)
    
    else:
        if anc == 'EUR':
            j._machine_type = 'n1-highmem-16'
        else:
            j.memory(str(memsize) + 'G')
        j.storage(str(storage) + 'G')
        
        genotype_name = get_geno_split_names([ancestry], dictout=True, plinkprefix=True)[ancestry]
        covariate_name = get_covar_split_names([ancestry], dictout=True, sex_specific=sex_specific)[ancestry]
        annot_name = get_annot_split_names([ancestry], dictout=True, n_annot=nbins)[ancestry]

        covarpath = b.read_input(path_covar + covariate_name)
        annotpath = b.read_input(path_annot + annot_name)

        if use_fuse:
            genopath = f'/{fusebucket}/genos/{genotype_name}'
            j.gcsfuse(fusebucket, '/'+fusebucket, read_only=True)
        else:
            genopath = b.read_input_group(**{x: path_geno + genotype_name + '.' + x 
                                            for x in ['bed', 'bim', 'fam']})
        if True: # TODO update this to include input flags
            command_rhe = f"""
                /RHE-mc/build/RHEmc \
                -g {genopath} \
                -c covar.tsv \
                -p pheno.tsv \
                -annot annot.txt \
                -k {random_vectors} \
                -jn {jackknife_blocks} \
                -o {j.h2out}
            """
        else:
            # single component RHE reg either runs out of memory or
            # doesn't work. Gets floating point error:
            # 0x00005571c3c03e86 in genotype::get_observed_pj(unsigned char const*) ()
            command_rhe = f"""
                /RHE-reg/build/RHE_reg \
                -g {genopath} \
                -c covar.tsv \
                -p pheno.tsv \
                -b {random_vectors}
            """
        command = f"""
            cp {phenotype_file} pheno.bgz
            gunzip -c pheno.bgz > pheno.tsv
            cp {annotpath} annot.bgz
            gunzip -c {annotpath} > annot.txt
            cat {covarpath} > covar.tsv

            {command_rhe}

            python3 {parser} \
                --file {j.h2out} \
                --pheno '{phenotype_id}{iter_suffix}' \
                --ancestry {ancestry} \
                --out {j.parsed_out}
        """
        j.command(command)
        b.write_output(j.h2out, path_results + 'raw/' + log_out)
    
    return b, j, j.parsed_out


def run_ancestry_sink(b, phenotype_id, ancestry_jobs, concatter):
    """ Runs ancestry sink, producing tables for final concatenation.
    """
    j = b.new_job(name='ancestry_sink_' + phenotype_id)
    j.image(IMAGE)
    tables = ','.join([res_file for _, (_, res_file) in ancestry_jobs.items()])
    command = f"""
        python {concatter} --tables {tables} --out {j.tab_out}
    """
    j.command(command)
    
    return j


def run_final_sink(b, ancestry_sinks, concatter, nlen, suffix='', output_file='tab_out'):
    """ Runs final sink to collect and concatenate all results.
    Implements interim sinks of size nlen, and then has one final sink.
    This is to workaround the issue with the submitted script being too long
    if we allow the final sink size to become unbounded.
    """
    interim_sinks = []
    get_interim_id = lambda idx: math.floor(idx/nlen)
    for this_id in range(0, get_interim_id(len(anc_sinks))+1):
        jobs_in_sink = [v for idx, v in enumerate(ancestry_sinks) if get_interim_id(idx) == this_id]
        interim_sink = b.new_job('interim_sink_' + str(this_id))
        interim_sink.image(IMAGE)
        interim_tables = ",".join([j.tab_out for j in jobs_in_sink])
        command_interim = f"""
            python {concatter} --tables {interim_tables} --out {interim_sink.tab_out}
            """
        interim_sink.command(command_interim)
        interim_sinks.append(interim_sink)

    if len(interim_sinks) == 1:
        final_sink = interim_sinks[0]
    else:
        final_sink = b.new_job('final_sink')
        final_sink.image(IMAGE)
        final_tables = ",".join([j.tab_out for j in interim_sinks])
        command_fin = f"""
            python {concatter} --tables {final_tables} --out {final_sink.tab_out}
            """
        final_sink.command(command_fin)
    bserv.write_output(final_sink.tab_out, f'{path_results}final_results{suffix}.tsv')

    return final_sink


def get_num_bins(args):
    """ Gets the number of bins that will be used in this analysis.
    """
    n_maf = 2 if args.maf_bins_2 else 5
    n_ld = args.num_ld_bins
    return n_maf * n_ld


def get_num_iter(args):
    """ Determines if iterations per phenotype must be performed (only enabled for testing).
    If so, formats iterations properly for use.
    """
    if args.n_iter is None:
        return [None]
    else: 
        return range(0, args.n_iter)


def _get_pheno_manifest_path_internal():
    return f'gs://ukb-diverse-pops/{loc}/phenotype_manifest.tsv.bgz'


def _import_manifest(use_tsv_manifest=True):
    if use_tsv_manifest:
        manifest = hl.import_table(_get_pheno_manifest_path_internal())
    else:
        manifest = hl.read_matrix_table(get_variant_results_path('full')).cols()
        annotate_dict = {}
        annotate_dict.update({'pops': hl.str(',').join(manifest.pheno_data.pop),
                              'num_pops': hl.len(manifest.pheno_data.pop)})
        for field in ['n_cases','n_controls','heritability']:
            for pop in ['AFR','AMR','CSA','EAS','EUR','MID']:
                new_field = field if field!='heritability' else 'saige_heritability' # new field name (only applicable to saige heritability)
                idx = manifest.pheno_data.pop.index(pop)
                field_expr = manifest.pheno_data[field]
                annotate_dict.update({f'{new_field}_{pop}': hl.if_else(hl.is_nan(idx),
                                                                       hl.missing(field_expr[0].dtype),
                                                                       field_expr[idx])})
        manifest = manifest.annotate(**annotate_dict)
        manifest = manifest.drop(manifest.pheno_data)
    if 'phenotype_id' not in list(manifest.row):
        manifest = manifest.annotate(phenotype_id = construct_phenotype_id(manifest))
    manifest = manifest.annotate(pheno_file = get_pheno_filename(manifest.phenotype_id),
                                 saige_file = convert_pheno_id_to_potential_saige(manifest.phenotype_id),
                                 pop_split = manifest.pops.split(','))
    return manifest.cache()


def _make_phenotype_dict(manifest, ancestries, n_include=None, random_phenos=None, suffix='', specific_pheno=None):
    """ Internal function to construct a phenotype dictionary. Used to iterate
    through phenotypes during pipeline construction.

    Parameters
    ----------
    manifest: `Hail Table`
    Imported manifest table. This function relies on modifications made to this table
    in `_import_manifest`.

    ancestries: `list`

    n_include: `int` or None
    This is primarily used for testing. If a non-None value is provided,
    this function will trim the dictionary to include n_include items from each
    trait type.

    random_phenos: `int` or None
    If not none, then will queue up a bunch of random phenotype jobs. Each of these will
    be run on all phenotypes.

    suffix: `str`
    Appended to the end of the results file to search for.

    specific_pheno: `list` or None
    A list of formatted phenotypes to subset the manifest to. If any are not found, an error will be thrown.
    """
    if random_phenos is not None:
        dct_out = {}
        for i in range(1, random_phenos+1):
            pheno_id = RANDOMPREF(i)
            pheno_file = get_pheno_filename(pheno_id)
            saige_file = convert_pheno_id_to_potential_saige(pheno_id)
            dct_out.update({pheno_id: (pheno_file, saige_file, ancestries)})
        return dct_out
    else:
        phenotype_id = manifest.phenotype_id.collect()
        pheno_file = manifest.pheno_file.collect()
        saige_file = manifest.saige_file.collect()
        anc_vec = [[anc for anc in anclist if anc in ancestries] 
                for anclist in manifest.pop_split.collect()]

        zipped_vals = list(zip(phenotype_id, pheno_file, saige_file, anc_vec))
        pheno_to_len = {phen: len(anclist) for phen, _, _, anclist in zipped_vals}
        dct_out = {id: (file, saigefile, anclist) for id, file, saigefile, anclist in zipped_vals if len(anclist) > 0}

        if n_include is not None:
            trait_types = list(manifest.aggregate(hl.agg.collect_as_set(manifest.trait_type)))
            list_of_pheno_lists = []
            for trait_type in trait_types:
                corresponding_phenos = manifest.filter(manifest.trait_type == trait_type).phenotype_id.collect()
                corr_phenos_slim = [pheno for pheno in corresponding_phenos if pheno_to_len[pheno] > 0]
                list_of_pheno_lists.append(corr_phenos_slim[0:min(n_include, len(corr_phenos_slim))])
            traits_to_keep = [item for pheno_list in list_of_pheno_lists for item in pheno_list]
            return {key: dct_out[key] for key in traits_to_keep}
        elif specific_pheno is not None:
            if not all([x in dct_out.keys() for x in specific_pheno]):
                missing_pheno = ', '.join([x for x in specific_pheno if x not in dct_out.keys()])
                raise ValueError('ERROR: the following phenotypes were not found: ' + missing_pheno)
            return {key: dct_out[key] for key in specific_pheno}
        else:
            return dct_out


def _read_pheno_rowtable():
    """ Reads the row Hail Table from the FULL phenotype MatrixTable.
    """
    return _read_pheno_data(False,load_full=True).rows()


def _read_pheno_data(checkpoint, load_full=False):
    """ Reads phenotype data, extending the stock phenotype matrix table
    by using a single phenotype ID and keeping only a single row, entry, and column field
    unless load_full is set to True.

    Returns
    -------
    `MatrixTable`
    """
    checkpoint_phenos = MT_TEMP_BUCKET + 'filtered_phenotype_data_temporary.mt'
    if load_full:
        return hl.read_matrix_table(get_ukb_pheno_mt_path())
    else:
        if hl.hadoop_exists(checkpoint_phenos):
            pheno_mt_string = hl.read_matrix_table(checkpoint_phenos)
        else:
            pheno_mt = hl.read_matrix_table(get_ukb_pheno_mt_path())
            pheno_mt_string = pheno_mt.annotate_cols(phenotype_id = construct_phenotype_id(pheno_mt)
                                     ).key_cols_by('phenotype_id')
            pheno_mt_string = pheno_mt_string.select_entries('both_sexes'
                                            ).select_cols().select_rows(
                                            ).rename({'both_sexes':'value'})
            if checkpoint:
                pheno_mt_string = pheno_mt_string.checkpoint(checkpoint_phenos)
        return pheno_mt_string


def _write_flat_without_key(ht, destination, delimiter, header):
    """ Helper function to write the contents of a Hail table to delimited flat file.
    This function extends the usual Hail Table export function, however guarentees
    preserved order upon unkeying and then writes the corresponding file.
    """
    ht2 = ht.order_by(*ht.key).drop(*ht.key)
    ht2.export(output=destination, header=header, delimiter=delimiter)


def compatiblify_phenotype_id(phenotype_id):
    """RHEmc throws errors on reading weird phenotype
    headers. This function resolves these characters.
    """
    to_replace = [' ', '|', '>', '<', '/', '\\']
    for i in to_replace:
        phenotype_id = phenotype_id.replace(i, '_')
    return phenotype_id


def _verify_args(args):
    """ Ensures that arguments are properly inputted.
    """
    if args.initialize_only and args.phenotype_only:
        raise argparse.ArgumentError(message='Only one of --initialize-only and --phenotype-only can be specified.')
    if args.n_only is not None and args.n_only <= 0:
        raise argparse.ArgumentError(message='If --n-only is set, it must be greater than 0.')
    if (args.initialize_only or args.phenotype_only) and args.read_previous_rhemc:
        raise argparse.ArgumentError(message='Cannot use --initialize-only or --phenotype-only with --read-previous-rhemc.')
    if (args.initialize_only or args.phenotype_only) and args.n_iter is not None:
        raise argparse.ArgumentError(message='Cannot use --initialize-only or --phenotype-only with --n-iter.')
    if args.random_phenotypes is not None and (args.initialize_only or args.read_previous_rhemc or (args.n_only is not None)):
        raise argparse.ArgumentError(message='Since --random-phenotypes is enabled, pipeline cannot be run with ' + 
                                     '--read-previous-rhemc, --n-only, or --initialize-only.')
    if args.specific_pheno is not None and (args.random_phenotypes is not None or args.initialize_only or (args.n_only is not None)):
        raise argparse.ArgumentError(message='Since --specific-pheno is enable, pipeline cannot be run with ' +
                                     '--n-only, --initialize-only, or --random-phenotypes.')


def _initialize_log(args):
    logging.info('RHEmc HE pipeline for UKB Pan Ancestry heritability analysis')
    logging.info('Rahul Gupta, 2021')
    logging.info(time.ctime())
    logging.info('-----------------------------------------')
    logging.info('Call:')
    logging.info('python rhemc_pipeline.py')
    for k,v in vars(args).items():
        logging.info(f'{"--" + k.replace("_", "-")} {v}')
    logging.info('-----------------------------------------')

    if args.initialize_only:
        logging.info('NOTE: Running with --initialize-only and thus only verifying covariate, annotation, and genotype files.')
    if args.n_only is not None:
        logging.info('NOTE: Will use at most ' + str(args.n_only) + ' traits from each category.')
    if args.phenotype_only:
        logging.info('NOTE: Running with --phenotype-only and thus will create all data files (including phenotype files).')
    if args.read_previous_rhemc:
        logging.info('NOTE: We will read previously completed RHEmc runs from ' + path_results + ' if available.')
    if args.random_phenotypes is not None:
        logging.info('NOTE: Running random phenotypes only. Will do ' + str(args.random_phenotypes) + ' phenos per ancestry.')
    if args.n_iter is not None:
        logging.info('NOTE: Running multiple duplicates of each penotype. Will do ' + str(args.n_iter) + ' runs for each phenotype-ancestry pair.')
    if args.specific_pheno is not None:
        logging.info('NOTE: Will only test the following phenotypes: ' + ', '.join(args.specific_pheno.split(',')))
    if args.sex_specific:
        logging.info('NOTE: Sex-specific flag enabled. Sex covariates will not be used.')
        if args.specific_pheno is None:
            logging.info('Because specific phenotypes were not specified for this sex-specific analysis, phenotypes will be filtered to those with pheno_sex != both_sexes.')


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--ancestries', default='AFR,AMR,CSA,EAS,EUR,MID',
                        help='Comma-delimited set of ancestries to include. Default is all 6.')
    parser.add_argument('--suffix', default='',
                        help='A suffix appended to all outputted results files (logs and the final file).')
    parser.add_argument('--verbose', action='store_true',
                        help='If enabled, will print status updates to std out.')
    parser.add_argument('--checkpoint', action='store_true',
                        help='If enabled, will checkpoint all individual level data in '+ \
                            'Hail format prior to outputting flat/binary files in check_indiv_files.')
    parser.add_argument('--read-previous-rhemc', action='store_true',
                        help='If enabled, will check the results folder for previous results from RHE-mc.' + \
                            ' If previous results are found, these will be pulled rather than running the full pipeline. ' +\
                            'Cannot be enabled if --initialize-only or --phenotype-only are specified.')
    parser.add_argument('--logging', default='rhemc_pipeline.log', 
                        help='Path for pipeline log file. Note that logs will be produced when ' + \
                            'creating phenotype files as well.')
    parser.add_argument('--use-fuse', action='store_true',
                        help='If enabled, will localize individual level data with GCS fuse.')                 
    parser.add_argument('--approx-quantiles', action='store_true',
                        help='If enabled, will use approximate quantiles for LD score bins.' + \
                            'Note that this is non-deterministic.')
    parser.add_argument('--num-ld-bins', default=5, type=int,
                        help='Number of LD bins to include. Created as quantiles.')
    parser.add_argument('--maf', default=0.01, type=float,
                        help='Minimum MAF cutoff for pipeline.')
    parser.add_argument('--maf-bins-2', action='store_true',
                        help='If enabled, will use two MAF bins. If disabled, will default to 5 MAF bins.')
    parser.add_argument('--sex-specific', action='store_true',
                        help='If enabled, will exclude sex covariates. This should be enabled only if the phenotypes to be run ' +
                        'are sex-stratified.')
    parser.add_argument('--use-tsv-manifest', action='store_true',
                        help='If enabled, will use a .tsv manifest. This is a legacy option, as using cols from the sumstats ' + \
                        'MatrixTable will be the most updated and is thus preferred.')

    parser.add_argument('--n-threads-pheno', default=4, type=int,
                        help='Number of threads per phenotype worker.')
    parser.add_argument('--mem-rhemc', default=8, type=int,
                        help='GB of memory allocated for each RHEmc worker.')
    parser.add_argument('--mem-rhemc-eur', default=104, type=int,
                        help='GB of memory allocated for each EUR RHEmc worker.')
    parser.add_argument('--store-rhemc', default=10, type=int,
                        help='GB of storage allocated for each non-EUR RHEmc worker. It is highly recommended to enable --use-fuse ' + \
                            'if the size of the genotype files are large, rather than expanding storage via this avenue.')
    parser.add_argument('--store-rhemc-eur', default=20, type=int,
                        help='GB of storage allocated for each EUR RHEmc worker. It is highly recommended to enable --use-fuse ' + \
                            'if the size of the genotype files are large, rather than expanding storage via this avenue.')
    parser.add_argument('--jackknife-blocks', default=100, type=int,
                        help='Number of jackknife blocks to use in RHEmc. 100 (default) or 22 recommended.' + \
                            ' Higher values result in more memory usage.')
    parser.add_argument('--random-vectors', default=10, type=int,
                        help='Number of random vectors to use for RHEmc. 10 (default) is recommended.')
    parser.add_argument('--random-phenotypes', type=int,
                        help='Construct and run RHEmc on a set of random phenotypes. If this flag is used, ' +
                            '--initialize-only, --read-previous-rhemc, and --n-only cannot be enabled. ' +
                            'Enabling this requires listing a number of random phenotypes to run. Random phenotypes ' +
                            ' are constructed from a Normal(0,1).')

    parser.add_argument('--initialize-only', action='store_true',
                        help='If enabled, this pipeline will only run check_indiv_files, which creates and checks for ' + \
                        'several important helper files including genotype MatrixTables, PLINK files, ' + \
                        '.fam file Hail Tables, the phenotype file MatrixTable, the covariate files, and ' + \
                        'annotation files. Note that this will not initialize phenotypes. To initialize files and create ' + \
                        'all phenotype files, enable --phenotype-only. This flag cannot be concurrently enabled with ' + \
                        '--phenotype-only.')
    parser.add_argument('--n-only', type=int,
                        help='Selects n of each type of trait for a test run.')
    parser.add_argument('--specific-pheno', 
                        help='Comma-delimited list of phenotypes to run. Error will be thrown if these phenotypes are not found' +\
                            ' in the manifest. See construct_phenotype_id() to see how phenotype IDs should be formatted. ' + \
                            'Cannot be enabled alongside --random-phenotypes, --n-only, or --initialize-only.')
    parser.add_argument('--n-iter', type=int,
                        help='To characterize run-to-run variability, specify a number of iterations here. If specified, an integer ' + \
                        'will be appended to phenotype id as "-iter_". Cannot be specified alongside --phenotype-only or --initialize-only.')
    parser.add_argument('--phenotype-only', action='store_true',
                        help='Only run phenotypes, without running RHEmc. Cannot be enabled if ' + \
                            '--initialize-only is enabled.')

    args = parser.parse_args()
    _verify_args(args)
    ancestries = parse_ancestries(args)
    nbins = get_num_bins(args)
    niter = get_num_iter(args)
    
    logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s", level='INFO', filename=args.logging)
    if args.verbose:
        logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))
    _initialize_log(args)
    
    # Ensures all files are present. Does not verify phenotype files; this is done in pheno workers
    check_indiv_files(ancestries, args, nbins)

    if not args.initialize_only and not args.phenotype_only and args.use_fuse:
        prepare_fuse_bucket(ancestries, nbins, args.sex_specific)

    if not args.initialize_only:
        # Log the files present in the phenos directory
        ls_dir = list_precomputed_pheno_files(ancestries)
        saige_dir = list_saige_pheno_files(ancestries)
        ls_previous_runs = list_completed_rhemc_logs()

        # Get manifest
        manifest = _import_manifest(args.use_tsv_manifest)
        if args.specific_pheno is None:
            if args.sex_specific:
                manifest = manifest.filter(manifest.pheno_sex != 'both_sexes')
            phenodict = _make_phenotype_dict(manifest, ancestries, args.n_only, 
                                             args.random_phenotypes, args.suffix,
                                             None)
        else:
            phenodict = _make_phenotype_dict(manifest, ancestries, args.n_only, 
                                            args.random_phenotypes, args.suffix,
                                            args.specific_pheno.split(','))

        # Create batch job
        backend = hb.ServiceBackend(billing_project='ukb_diverse_pops', bucket=bucket)
        bserv = hb.Batch(name="h2_rhemc_pan_ancestry", backend=backend)
        pheno_script = bserv.read_input(pheno_script_adr)
        cat_script = bserv.read_input(cat_script_adr)
        parsing_script = bserv.read_input(parsing_script_adr)
        
        # Run pipeline
        anc_sinks = []
        logging.info('Allocating pipeline. Analyzing ' + str(len(phenodict)) + ' phenotypes.')
        num_pheno_jobs = 0
        for pheno, (file, saigefile, anclist) in phenodict.items():
            # Figure out which ancestries do not have available phenotype files
            missing_ancestries = [anc for anc in anclist if file not in ls_dir[anc]]
            if (len(missing_ancestries) > 0) or (args.random_phenotypes is not None):
                pheno_job = True
                saige_miss = [anc for anc in anclist if saigefile not in saige_dir[anc]]
                bserv, j_pheno, map_dep = run_phenotype_job(b=bserv, phenotype_id=pheno, 
                                                            ancestries=anclist, 
                                                            use_saige=(len(saige_miss) == 0),
                                                            checkpoint=args.checkpoint, 
                                                            phenoscript=pheno_script,
                                                            n_threads=args.n_threads_pheno,
                                                            random=args.random_phenotypes is not None)
                num_pheno_jobs+=1
            else: 
                pheno_job = False # all phenotype files found!

            if not args.phenotype_only:
                ancestry_jobs = {}
                for anc in anclist:
                    if pheno_job:
                        this_pheno_file = map_dep[anc]
                    else:
                        this_pheno_file = bserv.read_input(path_pheno + anc + '/' + file)

                    mem_use = args.mem_rhemc_eur if anc == 'EUR' else args.mem_rhemc
                    storage_use = args.store_rhemc_eur if anc == 'EUR' else args.store_rhemc
                    
                    for niter_idx in niter:
                        iter_suff = construct_iter_suffix(niter_idx)
                        log_out_dir = get_result_log_split_names([anc], pheno, dictout=True, suffix=args.suffix)[anc] + iter_suff
                        read_previous = args.read_previous_rhemc and log_out_dir in ls_previous_runs
                        bserv, j_rhemc, rhemc_res = run_rhemc_job(bserv, pheno, anc, this_pheno_file,
                                                                read_previous=read_previous,
                                                                parser=parsing_script,
                                                                use_fuse=args.use_fuse,
                                                                memsize=mem_use,
                                                                storage=storage_use,
                                                                jackknife_blocks=args.jackknife_blocks,
                                                                random_vectors=args.random_vectors,
                                                                nbins=nbins,
                                                                suffix=args.suffix, iter=niter_idx,
                                                                sex_specific=args.sex_specific)
                        if niter is None:
                            anc_iter = anc
                        else:
                            anc_iter = anc + '_iter' + str(niter_idx)
                        
                        ancestry_jobs.update({anc_iter: (j_rhemc, rhemc_res)})
                
                anc_sink = run_ancestry_sink(bserv, pheno, ancestry_jobs, cat_script)
                anc_sinks.append(anc_sink) # these jobs collect ancestry level h2 estimates
        
        logging.info(str(num_pheno_jobs) + ' phenotype jobs successfully created. The rest ('+ \
                     str(len(phenodict)-num_pheno_jobs) + ') were previously computed.')
        if not args.phenotype_only:
            final_sink = run_final_sink(bserv, anc_sinks, cat_script, nlen=500, suffix=args.suffix)
            logging.info('RHEmc jobs and sinks successfully created.')
        
        logging.info('Starting pipeline...')
        bserv.run(verbose=False)
