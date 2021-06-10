#!/usr/bin/env python

__author__ = 'Rahul Gupta'

import re
import argparse
import pandas as pd

SCINOT = r'[0-9\-\+.enan]'

def internal_parse_args(parser):
    """ Parse incoming arguments, catch issues, and convert to lists.

    Keyword arguments:

    parser -- an ArgumentParser()

    Returns:
    dict of arguments: new values
    """
    args = parser.parse_args()
    if args.logs is None:
        raise ValueError('Must provide at least one log file.')
    if args.code is None:
        raise ValueError('Must provide phenotype code.')
    if args.ancestry_vec is None:
        raise ValueError('Must provide at least one ancestry.')
    if args.out is None:
        raise ValueError('Must provide output directory.')
    list_of_logs = args.logs.split(',')
    list_of_ancestry = args.ancestry_vec.split(',')
    if len(list_of_logs) != len(list_of_ancestry):
        raise ValueError('The list of logs and list of ancestries must be the same length.')
    args_dict = vars(args)
    args_dict.update({'logs': list_of_logs, 'ancestry_vec': list_of_ancestry})
    return args_dict


def import_log(log_directory):
    """ Import a log file from a provided directory.

    Keyword arguments:

    log_directory -- a path on disk.

    Returns:
    List of lines imported.
    """
    list_of_lines = []
    with open(log_directory, 'r') as reader:
        for line in reader:
            this_line = line.strip()
            list_of_lines.append(this_line)
    reader.close()
    return list_of_lines


def generate_table_row(log_file, ancestry, official_only, code):
    """ Takes an imported log and ancestry and converts it into a properly formatted pandas table.

    Keyword arguments:

    log_file -- output from import_log()

    ancestry -- a single ancestry code

    official_only -- a boolean indicating if all fields should be imported
    into the table, or only the official ones.

    Returns:
    dict of arguments: new values
    """
    # verify that ancestry is correct
    matches = [l for l in log_file if re.search('Searching for ancestry: ' + \
                                                ancestry, l)]
    if len(matches) == 0:
        raise ValueError('ALERT: Incorrect ancestry passed in for code ' + code +
                         '. Passed in value: ' + ancestry)
    dict_of_vals = {'ancestry': ancestry, 'phenotype_code': code}

    nrow_orig = num_cols = None
    for line in log_file:
        nrow_orig = _parse_single_term(nrow_orig, 'Original number of rows: ([0-9]*)',
                                       line, int)
        num_cols = _parse_single_term(num_cols, 'Found ([0-9]*) ancestry specific columns:',
                                      line, int)
    dict_of_vals.update({'original_nrow': nrow_orig, 
                         'ancestry_specific_ncols': num_cols})
    
    if dict_of_vals['ancestry_specific_ncols'] != 0:
        tf_boundary = [idx for idx, l in enumerate(log_file) if re.search('Now running LDSC in vanilla mode.',l)]
        log_file_official = log_file[(tf_boundary[0]+1):(len(log_file)+1)]
        log_file_unofficial = log_file[0:tf_boundary[0]]
        if not official_only:
            unofficial_dict = _parse_unofficial_log(log_file_unofficial)
            dict_of_vals.update(unofficial_dict)
        official_dict, error_str = _parse_official_log(log_file_official)
    else:
        if not official_only:
            unofficial_dict = _parse_unofficial_log(log_file)
            dict_of_vals.update(unofficial_dict)
        official_dict, _ = _parse_official_log(log_file)
        error_str = 'No ' + ancestry + '-specific columns found.'
    
    dict_of_vals.update(official_dict)
    if error_str is not None:
        dict_of_vals.update({'missing_data_note': error_str})
    return pd.DataFrame(dict_of_vals, index=[ancestry + ':' + code])


def _parse_single_term(ele, string, line, type_in, group=1):
    """ Stereotyped approach to quickly test for matches to a particular string.

    Keyword arguments:

    ele -- Currently saved element

    string -- The string to search for. Expects at least one group.

    line -- Current line being searched.

    type_in -- python type to coerce found results to.

    group -- the group number. Defaults to 1, but can be overridden.

    Returns:
    The extracted element of type determined by type_in, or np.nan
    """
    if ele is None:
        ele_s = re.search(string, line)
        if ele_s:
            return type_in(ele_s.group(group))
    else:
        return ele


def _parse_unofficial_log(list_lines):
    """ Takes an imported log and returns unofficial value dict.
    Returns:
    dict of arguments: new values
    """
    nrow_orig = num_cols = N = nrow_premunge = nrow_postmunge = None
    nrow_rem_na = nrow_rem_oob_p = nrow_rem_af = nrow_rem_low_conf = None
    mean_chi2 = lgc = max_chi2 = n_genome_sig = None

    for line in list_lines:
        if num_cols is not None and num_cols == 0:
            break
        N = _parse_single_term(N, 'INFO:root:Used N = ([0-9]*[.][0-9]*)[.]',
                               line, float)
        nrow_premunge = _parse_single_term(nrow_premunge, 'INFO:root:Original row count: ([0-9]*)',
                                           line, int)
        nrow_postmunge = _parse_single_term(nrow_postmunge, 'INFO:root:Final number of rows in file: ([0-9]*)',
                                            line, int)
        nrow_rem_na = _parse_single_term(nrow_rem_na, 'INFO:root:Rows removed due to NA in [A-Za-z,_ ]+: ([0-9]*)',
                                         line, int)
        nrow_rem_oob_p = _parse_single_term(nrow_rem_oob_p, 'INFO:root:Rows removed due to pval_.{3} out of bounds: ([0-9]*)',
                                            line, int)
        nrow_rem_af = _parse_single_term(nrow_rem_af, 'INFO:root:Rows removed due to af_(controls_)?.{3} below 0[.]01: ([0-9]*)',
                                         line, int, 2)
        nrow_rem_low_conf = _parse_single_term(nrow_rem_low_conf, 'INFO:root:Rows removed due to low confidence: ([0-9]*)',
                                               line, int)
        mean_chi2 = _parse_single_term(mean_chi2, 'INFO:root:Mean chi\\^2 = ([0-9.]*)', line, float)
        lgc = _parse_single_term(lgc, 'INFO:root:Lambda GC = ([0-9.]*)', line, float)
        max_chi2 = _parse_single_term(max_chi2, 'INFO:root:Max chi\\^2 = ([0-9.]*)', line, float)
        n_genome_sig = _parse_single_term(n_genome_sig, 'INFO:root:([0-9]*) genome-wide significant SNPs \\(some may have been removed by filtering\\)',
                                          line, int)
    return {'N': N, 'premunge_nrow': nrow_premunge, 'postmunge_nrow': nrow_postmunge,
            'nrow_with_NA_removed': nrow_rem_na, 'nrow_with_out_of_bound_p': nrow_rem_oob_p,
            'nrow_with_low_AF': nrow_rem_af, 'nrow_with_low_conf': nrow_rem_low_conf,
            'custom_mean_chi2': mean_chi2, 'custom_max_chi2': max_chi2,
            'custom_lambda_GC': lgc, 'num_genomewide_sig_postmunge': n_genome_sig}


def _parse_official_log(list_lines):
    """ Takes an imported log and returns official value dict.
    Returns:
    dict of arguments: new values
    """
    nrow_orig = nrow_ldscores = nrow_postmerge = None
    h2 = h2_sd = lgc = mean_chi2 = None
    intercept = intercept_sd = ratio = ratio_sd = miss_rat_less1 = miss_rat_gc = None

    for line in list_lines:
        nrow_orig = _parse_single_term(nrow_orig, 'Read summary statistics for ([0-9]*) SNPs',
                                       line, int)
        nrow_ldscores = _parse_single_term(nrow_ldscores, 'Read reference panel LD Scores for ([0-9]*) SNPs.',
                                      line, int)
        nrow_postmerge = _parse_single_term(nrow_postmerge, 'After merging with regression SNP LD, ([0-9]*) SNPs remain.',
                               line, int)
        h2 = _parse_single_term(h2, 'Total Observed scale h2: ('+ SCINOT+ '*) \\(',
                                           line, float)
        h2_sd = _parse_single_term(h2_sd, 'Total Observed scale h2: '+ SCINOT+ '* \\(('+ SCINOT+ '*)\\)',
                                           line, float)
        mean_chi2 = _parse_single_term(mean_chi2, 'Mean Chi\\^2: ('+ SCINOT+ '*)', line, float)
        lgc = _parse_single_term(lgc, 'Lambda GC: ('+ SCINOT+ '*)', line, float)
        intercept = _parse_single_term(intercept, 'Intercept: ('+ SCINOT+ '*) \\(',
                                           line, float)
        intercept_sd = _parse_single_term(intercept_sd, 'Intercept: '+ SCINOT+ '* \\(('+ SCINOT+ '*)\\)',
                                          line, float)
        ratio = _parse_single_term(ratio, 'Ratio: ('+ SCINOT+ '*) \\(',
                                   line, float)
        ratio_sd = _parse_single_term(ratio_sd, 'Ratio: '+ SCINOT+ '* \\(('+ SCINOT+ '*)\\)',
                                      line, float)
        miss_rat_less1 = _parse_single_term(miss_rat_less1, 'Ratio: (NA) \\(mean chi\\^2 < 1\\)',
                                            line, str)
        miss_rat_gc = _parse_single_term(miss_rat_gc, 'Ratio < (0) \\(usually indicates GC correction\\)',
                                         line, str)
    
    data_dict = {'ldsc_original_nrow': nrow_orig, 'ldsc_ldscores_nrow': nrow_ldscores,
                 'ldsc_postmerge_nrow': nrow_postmerge, 'ldsc_h2': h2, 'ldsc_h2_sd': h2_sd,
                 'ldsc_intercept': intercept, 'ldsc_intercept_sd': intercept_sd,
                 'ldsc_ratio': ratio, 'ldsc_ratio_sd': ratio_sd,
                 'ldsc_mean_chi2': mean_chi2, 'ldsc_lambda_GC': lgc}
    
    if miss_rat_less1 is not None:
        error_str = 'No ratio as mean chi^2 < 1'
    elif miss_rat_gc is not None:
        error_str = 'Ratio < 0, possibly indicating GC correction'
    else:
        error_str = None

    return data_dict, error_str


parser = argparse.ArgumentParser()
parser.add_argument('--logs', type=str, default=None,
                    help='Comma delimited list of logs to import. Must have same ' + \
                    'number of elements as the input to --ancestry-vec.')
parser.add_argument('--official-only', type=bool, default=False,
                    help='Flag to enable official fields only (eg those provided directly ' +\
                         'by ldsc itself.')
parser.add_argument('--code', type=str, default=None,
                    help='Single value representing the phenotype code analyzed.')
parser.add_argument('--ancestry-vec', type=str, default=None,
                    help='Comma delimited list of ancestries. Must have the same ' + \
                         'number of elements as the input to --logs.')
parser.add_argument('--out', type=str, default=None,
                    help='Output directory.')

if __name__ == '__main__':
    args = internal_parse_args(parser)
    df_list = []
    for log, anc in zip(args['logs'], args['ancestry_vec']):
        log_file = import_log(log)
        log_tab = generate_table_row(log_file, anc, 
                                     args['official_only'], args['code'])
        df_list.append(log_tab)
    final_table = pd.concat(df_list)
    final_table.to_csv(sep='\t', path_or_buf=args['out'], index=False)
