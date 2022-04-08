__author__ = 'Rahul Gupta'

import hailtop.batch as hb
import hail as hl
hl.init(spark_conf={'spark.hadoop.fs.gs.requester.pays.mode': 'CUSTOM',
                   'spark.hadoop.fs.gs.requester.pays.buckets': 'ukb-diverse-pops-public',
                   'spark.hadoop.fs.gs.requester.pays.project.id': 'ukbb-diversepops-neale'})
from hail import hadoop_exists, hadoop_ls
import argparse
import os, re, math, sys
import pandas as pd
import numpy as np

from ukbb_pan_ancestry.heritability.ldsc_ukbb_div_pops_constants import *
from ukbb_pan_ancestry.resources.results import get_variant_results_path
from ukbb_pan_ancestry.export_results import get_pheno_id

# Set up locations for python scripts
munging_script_adr = f'gs://{bucket}/munge_manual.py'
parsing_script_adr = f'gs://{bucket}/parse_ldsc_log.py'
cat_script_adr = f'gs://{bucket}/concat_tables.py'
check_script_adr = f'gs://{bucket}/check_munging_complete.py'
IMAGE = 'gcr.io/ukbb-diversepops-neale/rgupta_ldsc'
PROJECT = 'ukb_diverse_pops'
REQPAYS = 'ukbb-diversepops-neale'


def get_N(manifest_row, ancestry):
    """ Get sample size from a particular row of the manifest.

    Parameters
    ----------
    manifest_row : :obj:`DataFrame`

    Returns
    -------
    :obj: `float`
    """
    nca = float(manifest_row['n_cases_'+ancestry])
    nco = float(manifest_row['n_controls_'+ancestry])
    return np.nansum([nca, nco])


def read_ld_score(b, ancestries, anc_to_ldscore, expect_bgz, prefix):
    """ Read LD score files. Expects that there is one file per ancestry, *not* split per chromosome.

    Parameters
    ----------
    b : :obj:`Batch`
    ancestries : `list` 
        Ancestry abbreviations compatible with `anc_to_ldscore`
    anc_to_ldscore : `function`
        Function that outputs an ldscore address taking in an ancestry abbreviation
    expect_bgz : `bool`
        If True, will expect bgz suffix for .ldscore files. Otherwise will expect .gz.
    prefix : `str`
        Prefix of file names that correspond to all ldscore files. For example, for
        the file ld.ldscore.gz, the prefix is ld.

    Returns
    -------
    :obj: `Batch`, :obj: `Dictionary`mapping ancestries to ld score files
    """
    ldscore_file_dict = {}
    ldsc_suffix = f'.ldscore.{"b" if expect_bgz else ""}gz'
    for anc in ancestries:
        this_adr = anc_to_ldscore(anc)
        ig = b.read_input_group(**{'l2.ldscore.gz': this_adr + '.' + prefix + ldsc_suffix,
                                   'l2.M_5_50': this_adr + '.' + prefix + '.M_5_50',
                                   'l2.M': this_adr + '.' + prefix + '.M'})
        ldscore_file_dict.update({anc: ig})
    return b, ldscore_file_dict


def preallocate_munged_stats(bucket, ancestries):
    """ To speed up the submission process in which gs:// is checked for munged summary statistcs,
    we implement a preallocation data structure such that the cloud only needs to be queried once.
    Assumes that data is stored in bucket/ancestry_code/munged_sumstats/.

    Parameters
    ----------
    bucket : :obj:`str`
        Name of the bucket in which munged summary statistics were outputted.
    ancestries : :obj:`list`
        List of ancestry abbreviations to search.

    Returns
    -------
    :obj: `DataFrame` with the manifest, :obj: `list` containing the prefix types (indicating phenotype types)
    """
    dir_tail = [os.path.basename(item['path']) for item in hadoop_ls(bucket) if item['is_dir']]
    ancestries_keep = [anc for anc in ancestries if anc in dir_tail]
    return {anc: [item['path'] for item in hadoop_ls(bucket + anc + '/munged_sumstats/')] for anc in ancestries_keep}


def pull_pheno_manifest(main_dir, flat_file_location, specific_pheno, tsv_manifest):
    """ Obtain phenotype manifest from local source. This function will also produce 
    a list of prefixes corresponding to phenotype types. Further, this function pre-processes
    the table to include a column for file, which is the filename without its address. Will
    only include phenotypes that are found at the inputted location for flat files.

    Parameters
    ----------
    main_dir : :obj:`str`
        Local directory which contains the `phenotype_manifest.tsv`.

    tsv_manifest: :obj:`bool`
        If true, will use the .tsv file version of the manifest. Ohterwise will obtain manifest from columns of sumstat matrixtable.

    Returns
    -------
    :obj: `DataFrame` with the manifest, :obj: `list` containing the prefix types (indicating phenotype types)
    """
    available_files = [os.path.basename(loc['path']) for loc in hadoop_ls(flat_file_location) if not loc['is_dir']]
    
    if tsv_manifest:
        pheno_manifest = pd.read_csv(main_dir + 'phenotype_manifest.tsv', "\t")
    else:
        pheno_manifest = hl.read_matrix_table(get_variant_results_path('full')).cols()
        annotate_dict = {}
        annotate_dict.update({'pops': hl.str(',').join(pheno_manifest.pheno_data.pop),
                              'num_pops': hl.len(pheno_manifest.pheno_data.pop)})
        for field in ['n_cases','n_controls','heritability']:
            for pop in ['AFR','AMR','CSA','EAS','EUR','MID']:
                new_field = field if field!='heritability' else 'saige_heritability' # new field name (only applicable to saige heritability)
                idx = pheno_manifest.pheno_data.pop.index(pop)
                field_expr = pheno_manifest.pheno_data[field]
                annotate_dict.update({f'{new_field}_{pop}': hl.if_else(hl.is_nan(idx),
                                                                       hl.missing(field_expr[0].dtype),
                                                                       field_expr[idx])})
        annotate_dict.update({'filename': get_pheno_id(tb=pheno_manifest)+'.tsv.bgz'})
        pheno_manifest = pheno_manifest.annotate(**annotate_dict)
        pheno_manifest = pheno_manifest.drop(pheno_manifest.pheno_data)
        pheno_manifest = pheno_manifest.to_pandas()

    if specific_pheno is not None:
        key_cols = ['trait_type', 'phenocode', 'pheno_sex', 'coding', 'modifier']
        specific_file = pd.read_csv(specific_pheno, '\t', dtype={x:str for x in key_cols})
        specific_file = specific_file.fillna('')
        if all([x in list(specific_file.columns) for x in key_cols]):
            specific_file = specific_file[key_cols]
            pheno_manifest = specific_file.merge(pheno_manifest, on = key_cols, how='left')
        else:
            raise ValueError('Argument to --specific-pheno must contain key columns.')

    pheno_manifest.loc[:, "file"] = pheno_manifest["filename"].apply(lambda x: os.path.basename(x))
    pheno_manifest_f = pheno_manifest.loc[pheno_manifest['file'].isin(available_files)]
    unique_pref = set([this_name.split("-")[0] for this_name in pheno_manifest_f['file']])
    return pheno_manifest_f, unique_pref


def check_mung_exists(anc, code, output_bucket, prealloc=None):
    """ Checks if both products of correct munging exists. Does not test if the mugning was successful, rather
    just assesses if the munging was completed.

    Parameters
    ----------
    anc : :obj:`str`
        Ancestry prefix.
    code : :obj:`str`
        Phenotype code.
    output_bucket : :obj:`str`
        The bucket in which the mugning output is stored.

    Returns
    -------
    :obj: `bool`
    """
    file_loc, log_loc = get_mung_locs(anc, code, output_bucket)

    if prealloc is None:
        mung_file_tf = hadoop_exists(file_loc)
        mung_log_tf = hadoop_exists(log_loc)
    else:
        # if prealloc is provided, use this much faster method
        found_ancestries = list(prealloc.keys())
        if anc in found_ancestries:
            mung_file_tf = file_loc in prealloc[anc]
            mung_log_tf = log_loc in prealloc[anc]
        else:
            mung_file_tf = False
            mung_log_tf = False

    # if either the log or the file are found, returns false
    return mung_file_tf and mung_log_tf


def get_mung_locs(anc, code, output_bucket):
    """ Convenience function to obtain the expected locations for munged scripts.

    Parameters
    ----------
    anc : :obj:`str`
        Ancestry prefix.
    code : :obj:`str`
        Phenotype code.
    output_bucket : :obj:`str`
        The bucket in which the mugning output is stored.

    Returns
    -------
    :obj: `str` with file location, :obj:`str` with log location
    """
    file_loc = output_bucket + anc + '/munged_sumstats/' + code + '.sumstats.gz'
    log_loc = output_bucket + anc + '/munged_sumstats/' + code + '.log'
    return file_loc, log_loc


def pipeline_setup_env(parsing_script):
    """ Set up munging/ldsc job. Install the proper LDSC, activate conda environment, and 
    ensure that the parsing script is named properly so that other scripts dependant on it 
    can still load it.

    Parameters
    ----------
    parsing_script : :obj:`InputResourceFile`
        Location of the parsing script as produced by the Batch object.

    Returns
    -------
    :obj: `str` with command to be run for setting up the job.
    """
    command_setup = f"""
        source activate ldsc
        cd ldsc
        git checkout 3d0c446
        cd ..
        cp {parsing_script} {os.path.dirname(parsing_script) + 'parse_ldsc_log.py'}"""
    return command_setup


def munge_sumstats(j, ancestry, ld_scores, sumstat_file, N, munging_script):
    """ Set up commands for munging summary statsitcs. This command assumes that summary statistics
    files are of the pan-ancestry variety and as such will subset the summary statistics to only
    look at the relevant ancestry groups. Will also generate a log to track munging steps while also
    keeping only variants that are found in the LD score file. Several additional munging steps are taken
    as part of `munge_manual.py` including:

    - Merge CHR/POS/REF/ALT into a single column (SNP)
    - Ensure ref is A1
    - Add a sample size column (N)
    - Produce Z scores from p-values signed by betas.
    - Remove any rows with NAs in pvalue, AF, beta, low_confidence columns
    - Remove low confidence variants
    - Remove variants with AF < 0.01
    - Remove variants with out of bounds p-values

    Will return the following line if ancestry specific columns are not found:
    `5 (or 6 in the case of case/control) ancestry specific columns not found; exiting.`

    Parameters
    ----------
    j : :obj: `Job`

    ancestry: :obj: `str`
        Ancestry abbreviation that will be used to subset the summary statistcs.
    ld_scores: :obj: `ResourceGroup`
        A resource group with the `.ldscore.gz`, `.M`, and `.M.5.50` files required for
        sumstat munging.
    sumstat_file: :obj: `InputResourceFile`
        Location of the summary statistics file, with format `.tsv.gz`.
    N : :obj:`float`
        Sample size.
    munging_script : :obj: `InputResourceFile`
        Location of the munging script, `munge_manual.py`

    Returns
    -------
    Tuple:

    :obj: `str`
        Command to be run for performing munging.
    :obj: `Job`

    :obj: `str`
        Location of the sumstats file to be outputted.
    :obj: `ResourceFile`
        Location of the munging log to be outputted.
    """
    j = j.storage('20Gi')
    j = j.memory('8Gi')
    command_pre_script = f"""
        gunzip -c {ld_scores}.l2.ldscore.gz > {j.extracted_ldscore}
        gunzip -c {sumstat_file} > {j.extracted_sumstats}
        touch {j.logout}

        awk 'FNR==NR {{ a[$2]; next }} {{ b=$1":"$2":"$3":"$4; if (b in a) {{ print }} }}' {j.extracted_ldscore} {j.extracted_sumstats} > {j.intermediate_sumstats_1}
        cat <(head -n1 {j.extracted_sumstats}) {j.intermediate_sumstats_1} > {j.intermediate_sumstats_2}
        printf "Extracted sumstats and trimmed by variants present in LD score file.\n" >> {j.logout}
        printf "Original number of rows: %s\n" $(wc -l < "{j.extracted_sumstats}") >> {j.logout}
        printf "New number of rows after filtering by LD score file: %s\n" $(wc -l < "{j.intermediate_sumstats_2}") >> {j.logout}

        head -n1 {j.intermediate_sumstats_2} | tr "\t" "\n" | grep {ancestry} | cat > {j.headers_ancestry_specific}
        head -n1 {j.intermediate_sumstats_2} | tr "\t" "\n" | head -n4 > {j.headers_all}

        printf "Searching for ancestry: %s\n" {ancestry} >> {j.logout}
        export nlines=$(wc -l < "{j.headers_ancestry_specific}")
        printf "Found %s ancestry specific columns:\n" "${{nlines}}" >> {j.logout}
        cat {j.headers_ancestry_specific} >> {j.logout}"""
    command_condition = f"""
        if [ "${{nlines}}" -eq 5 ] || [ "${{nlines}}" -eq 6 ]; then
            cat {j.headers_all} {j.headers_ancestry_specific} > {j.headers_final}
            cut -f $(paste -sd, <(cut -f1 -d: <(grep -Fxn "$(<{j.headers_final})" < <(head -n1 {j.intermediate_sumstats_2} | tr "\t" "\n")))) {j.intermediate_sumstats_2} > {j.intermediate_sumstats_3}
            printf "Trimmed sumstats files to contain only %s columns.\n" {ancestry} >> {j.logout}
            python {munging_script} --sumstats {j.intermediate_sumstats_3} --N {N} --out {j.final_sumstats_file} --logfile {j.log_tmp}
            cp {j.final_sumstats_file} {j.final_sumstats_file + ".sumstats.gz"}
            printf -- "----------------------------------" >> {j.logout}
            printf "\n" >> {j.logout}
            cat {j.log_tmp} >> {j.logout}
        else
            printf "5 (or 6 in the case of case/control) ancestry specific columns not found; exiting.\n" >> {j.logout}
            touch {j.final_sumstats_file}
        fi"""
    sumstat_out = j.final_sumstats_file
    command_out = command_pre_script+command_condition
    return command_out, j, sumstat_out, j.logout


def run_ancestry_specific_job(b, ancestry, code, address, ld_scores, ld_weights, 
                              suffix,
                              parsing_script, munging_script,
                              check_script, prealloc=None, stratified=False, rem_maxchisq=False):
    """ Quarterbacks the ancestry specific job. This entails munging/obtaining the munged summary
    statistics for a particular ancestry and phenotype code and running ldsc. We need the parsing
    script because one of its functions is used by the check script.

    Parameters
    ----------
    b : :obj: `Batch`

    ancestry: :obj: `str`
        Ancestry abbreviation that will be used to subset the summary statistcs.
    code: :obj: `str`
        Phenotype code to test.
    address: :obj: `str`
        The address of the summary statistics to obtain.
    parsing_script : :obj: `InputResourceFile`
        Location of the parsing script script, `parse_ldsc_log.py`
    munging_script : :obj: `InputResourceFile`
        Location of the munging script, `munge_manual.py`
    check_script : :obj: `InputResourceFile`
        Location of the check script, `check_munging_complete.py`
    stratified : `bool`
        Determines if stratified LDSC is to be run.
    rem_maxchisq : `bool`
        IF enabled, will set maxchisq to 9999 in order to eliminate this cutoff.

    Returns
    -------
    Tuple:

    :obj: `Batch`

    :obj: `Job`

    """
    j = b.new_job(name=ancestry + '_' + code)
    j.declare_resource_group(outfiles={"log": '{root}.log'})
    j.image(IMAGE)
    command = pipeline_setup_env(parsing_script)
    stat_loc, log_loc = get_mung_locs(ancestry, code, output_bucket)
    already_munged = check_mung_exists(ancestry, code, output_bucket, prealloc=prealloc)
    if already_munged:
        # if the munging has been performed:
        local_stats = b.read_input_group(**{'sumstats.gz': stat_loc})['sumstats.gz']
        local_log = b.read_input_group(**{'log': log_loc})['log']
        local_stats_use = local_stats
    else:
        # if the munging has not been performed:
        sumstat_file = bserv.read_input_group(**{'tsv.bgz': address})['tsv.bgz']
        command_app, j, local_stats, local_log = munge_sumstats(j, ancestry, ld_scores, sumstat_file, this_N, munging_script)
        command += command_app
        local_stats_use = local_stats + ".sumstats.gz"
    
    chisqsuff = '--chisq-max 9999' if rem_maxchisq else ''
    # now that the munged stats have been obtained, run LDSC
    print_coef = ' --print-coefficients' if stratified else ''
    mode = 'stratified' if stratified else 'vanilla'
    command += f"""
        cp {local_log} {j.logout_with_ldsc}
        if [ $(python {check_script} --log {local_log}) == "Complete" ]; then
            printf "\nNow running LDSC in {mode} mode.\n" >> {j.logout_with_ldsc}
            python /ldsc/ldsc.py --h2 {local_stats_use} --ref-ld {ld_scores} --w-ld {ld_weights} --out {j.h2_log} {print_coef} {chisqsuff}
            cat {j.h2_log}.log >> {j.logout_with_ldsc}
        fi
        """
    j.command(command)

    # save outputs to disk
    if not already_munged:
        bserv.write_output(local_log, log_loc)
        bserv.write_output(local_stats, stat_loc)
    bserv.write_output(j.logout_with_ldsc, f'{output_bucket}{ancestry}/{code}{suffix}.log')
    return b, j


def run_ancestry_sink(b, ancestries_used, ancestry_jobs, code, parsing_script):
    """ Generates a sink job to collate logs from ancestry specific jobs. As such, will
    output results for a single trait across all possible ancestries.

    Parameters
    ----------
    b : :obj: `Batch`

    ancestries_used: :obj: `list`
        List of ancestry abbreviations that were analyzed.
    ancestry_jobs: :obj: `list`
        List of ancestry specific jobs generated by `run_ancestry_specific_job`
    code: :obj: `str`
        Phenotype code.
    parsing_script : :obj: `InputResourceFile`
        Location of the parsing script, `parse_ldsc_log.py`

    Returns
    -------
    :obj: `Job`

    """
    anc_sink = b.new_job('ancestry_sink_' + code.replace("'", '_'))
    anc_sink.image(IMAGE)
    anc_log_files = ",".join([j.logout_with_ldsc for j in ancestry_jobs])
    quo_log = '"' + anc_log_files + '"' if re.search("'",anc_log_files) else anc_log_files
    quo_code = '"' + code + '"' if re.search("'",code) else "'" + code + "'"
    command_anc_sink = f"""
        source activate ldsc
        python {parsing_script} --logs {quo_log} --code {quo_code} --ancestry-vec {",".join(ancestries_used)} --out {anc_sink.tab_out}
        """
    anc_sink.command(command_anc_sink)
    return anc_sink


def parse_ancestries(args):
    return args.ancestries.split(',')


def create_map_anc_ldscore(args, anc_to_ldscore):
    if args.stratified is None:
        return anc_to_ldscore, False
    else:
        return lambda anc: args.stratified.format(anc), True


def create_anc_ldscore_vanilla(args):
    return lambda anc: args.vanilla_ld.format(anc)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--ancestries', default='AFR,AMR,CSA,EAS,EUR,MID', type=str,
                        help='Comma-delimited set of ancestries to include. Default is all 6.')
    parser.add_argument('--suffix', default='', type=str,
                        help='A suffix appended to all outputted results files (logs and the final file).')
    parser.add_argument('--use-tsv-manifest', action='store_true',
                        help='If enabled, will use a .tsv manifest. This is a legacy option, as using cols from the sumstats ' + \
                        'MatrixTable will be the most updated and is thus preferred.')

    parser.add_argument('--vanilla-ld', default='gs://rgupta-ldsc/ld/UKBB.{}', type=str,
                        help='This argument points to the original complete in-sample LD scores. ' + \
                            'If stratified is disabled, these will be used for LDSC. Otherwise, for ' + \
                            'S-LDSC these will be used as weights.')
    parser.add_argument('--stratified', type=str,
                        help='If enabled, will use stratified LD scores for analysis. ' + \
                        'LD scores will be obtained from the provided gs:// address. Assumes that ' + \
                        'files are pooled across all chromosomes. The format must be ' + \
                        'similar to those requested for LDSC, namely that the full gs address is provided ' + \
                        'excluding the file extension suffix and {} is used for any location in the address ' + \
                        'where ancestry code should be substituted. For example, provide "UKBB.{}" for the file family ' + \
                        '"UKBB.AFR.l2.ldscore.gz" and "UKBB.AFR.l2.M".')
    parser.add_argument('--remove-maxchisq', action='store_true',
                        help='If enabled, will eliminate max chisq filter by setting --chisq-max = 9999.')

    parser.add_argument('--n-only', type=int,
                        help='Runs only the first n traits. Unlike the RHEmc pipeline, this ' + \
                        'does not attempt to run n traits from each trait type.')
    parser.add_argument('--specific-pheno', type=str,
                        help='Local path to specific phenotypes. Must contain the 5 columns: ' + \
                            'trait_type, phenocode, pheno_sex, coding, modifier. Will be used to filter the manifest internally. ')

    args = parser.parse_args()
    ancestries = parse_ancestries(args)
    anc_to_ldscore = create_anc_ldscore_vanilla(args)
    map_ancestry_ldscore, stratified = create_map_anc_ldscore(args, anc_to_ldscore)
    
    # Read in data paths ancestry list per phenotype
    pheno_manifest, unique_pref = pull_pheno_manifest(project_dir, flat_file_location, args.specific_pheno, args.use_tsv_manifest)

    # Set up Batch backend and input files/scripts
    backend = hb.ServiceBackend(billing_project=PROJECT, bucket=bucket)
    bserv = hb.Batch(name="h2_ldsc_pan_ancestry", backend=backend)
    munging_script = bserv.read_input(munging_script_adr)
    parsing_script = bserv.read_input(parsing_script_adr)
    cat_script = bserv.read_input(cat_script_adr)
    check_script = bserv.read_input(check_script_adr)

    # Get LD score files
    bserv, ldscore_file_dict = read_ld_score(b=bserv, ancestries=ancestries,
                                             anc_to_ldscore=map_ancestry_ldscore,
                                             expect_bgz=False, prefix='l2')
    bserv, ldscore_weight_dict = read_ld_score(b=bserv, ancestries=ancestries,
                                               anc_to_ldscore=anc_to_ldscore,
                                               expect_bgz=False, prefix='l2')

    # Preallocate for search of gs://
    prealloc = preallocate_munged_stats(output_bucket, ancestries)
    
    # trial_dir = '/Volumes/rahul/Projects/2020_ukb_diverse_pops/Experiments/200715_compare_h2_results_gnomad_ldsc/'
    # phenos_to_pull = list(pd.read_csv(trial_dir + 'Data/pheno_list_trial.txt', sep='\t', names=['phenos']).phenos)
    # phenos_to_use = pheno_manifest.loc[pheno_manifest.phenocode.map(lambda x: x in phenos_to_pull),]
    # address_list = [flat_file_location + fl for fl in phenos_to_use.file]
    # finish_pheno = ["prescriptions-acetylcholinesterase_inhibitor|Alzheimer's-both_sexes",
    #                 "prescriptions-adamantane|Parkinson's|anti-viral-both_sexes",
    #                 "prescriptions-dopamine_pro-drug___decarboxylase_inhibitor|Parkinson's-both_sexes"]
    # finish_pheno_f = [ph + '.tsv.bgz' for ph in finish_pheno]
    # address_list = [flat_file_location + fl for fl in pheno_manifest.file if fl in finish_pheno_f]

    address_list = [flat_file_location + fl for fl in pheno_manifest.file]
    address_list_iter = address_list if args.n_only is None else address_list[0:args.n_only]

    anc_sinks = []
    for address in address_list_iter:
        term = os.path.basename(address)
        code = term.replace('.tsv.bgz','') # updated to properly obtain unique IDs
        this_row = pheno_manifest.loc[pheno_manifest['file'] == term, :]
        if this_row.shape[0] == 0:
            raise ValueError('ERROR: File ' + term +
                            ' not found in phenotype manifest.')
        elif this_row.shape[0] > 1:
            raise ValueError('ERROR: File ' + term +
                            ' found multiple times in phenotype manifest.')
        populations = list(this_row.pops)[0].split(',')
        ancestries_used = []
        ancestry_jobs = []
        for ancestry in ancestries:
            if ancestry in populations:
                ancestries_used.append(ancestry)
                this_N = get_N(this_row, ancestry)
                ld_scores = ldscore_file_dict[ancestry]
                ld_weights = ldscore_weight_dict[ancestry]
                
                bserv, anc_j = run_ancestry_specific_job(bserv, ancestry, code, address,
                                                         ld_scores, ld_weights, args.suffix,
                                                         parsing_script, munging_script, 
                                                         check_script, prealloc=prealloc,
                                                         stratified=stratified, 
                                                         rem_maxchisq=args.remove_maxchisq)
                
                ancestry_jobs.append(anc_j)

        if len(ancestries_used) > 0:
            anc_sink = run_ancestry_sink(bserv, ancestries_used, ancestry_jobs, code, parsing_script)
            anc_sinks.append(anc_sink)
    
    # Implement interim sinks of size 500, and then have one final sink.
    # This is to workaround the issue with the submitted script being too long
    # if we allow the final sink size to become unbounded.
    interim_sinks = []
    get_interim_id = lambda idx: math.floor(idx/500)
    for this_id in range(0, get_interim_id(len(anc_sinks))):
        jobs_in_sink = [v for idx, v in enumerate(anc_sinks) if get_interim_id(idx) == this_id]
        interim_sink = bserv.new_job('interim_sink_' + str(this_id))
        interim_sink.image(IMAGE)
        interim_tables = ",".join([j.tab_out for j in jobs_in_sink])
        command_interim = f"""
            source activate ldsc
            python {cat_script} --tables {interim_tables} --out {interim_sink.tab_out}
            """
        interim_sink.command(command_interim)
        interim_sinks.append(interim_sink)

    final_sink = bserv.new_job('final_sink')
    final_sink.image(IMAGE)
    if len(interim_sinks) == 0:
        final_tables = ",".join([j.tab_out for j in anc_sinks])
    else:
        final_tables = ",".join([j.tab_out for j in interim_sinks])
    command_fin = f"""
        source activate ldsc
        python {cat_script} --tables {final_tables} --out {final_sink.tab_out}
        """
    final_sink.command(command_fin)
    bserv.write_output(final_sink.tab_out, f'{output_bucket}final_results{args.suffix}.tsv')

    bserv.run(verbose=False)
