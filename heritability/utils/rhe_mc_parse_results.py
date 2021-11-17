__author__ = 'Rahul Gupta'

import argparse
import re
import pandas as pd


SCINOT = r'[0-9\-\+.enan]' # regex to pick up an element of scientific notation


def read_file(loc):
    with open(loc) as text_file:
        list_of_lines_pre = text_file.readlines()
    list_of_lines = [re.sub(r'\n$', '', x) for x in list_of_lines_pre]

    return list_of_lines


def _get_two_matches_per_lines(list_of_lines, regex, transform=lambda x:x):
    match_pre = [re.search(regex, f) for f in list_of_lines 
                 if re.search(regex, f)]
    match_out = [(transform(x.group(1)), transform(x.group(2))) for x in match_pre]
    return match_out


def _get_first_match(list_of_lines, regex, transform=lambda x:x):
    match_pre = [re.search(regex, f) for f in list_of_lines 
                 if re.search(regex, f)][0].group(1)
    match_out = transform(match_pre)

    return match_out


def _get_first_idx(list_of_lines, regex):
    idxes = [idx for idx,x in enumerate(list_of_lines) if re.match(regex, x)]
    return idxes[0]


def count_inputs(list_of_lines, detailed=False):
    bins = [f for f in list_of_lines if re.search(' SNPs in [0-9]{1,2}-th bin', f)]
    nbins = len(bins) + 1
    num_snps = _get_first_match(list_of_lines, r'^SNP in bim file ([0-9]{1,})', transform=int)
    num_covariates = _get_first_match(list_of_lines, r'Read in ([0-9]{1,}) Covariates..', transform=int)
    num_indivs = _get_first_match(list_of_lines, r'Number of Indvs :([0-9]{1,})', transform=int)
    num_sel_snps = _get_first_match(list_of_lines, r'Number of selected SNPs w.r.t  annot file : ([0-9]{1,})', transform=int)
    num_snps_per_block = _get_first_match(list_of_lines, r'Number of SNPs per block : ([0-9]{1,})', transform=int)
    snps_per_bin = [_get_first_match(list_of_lines, r'^([0-9]{1,}) SNPs in ' + str(idx) + r'-th bin', transform=int) 
                    for idx in range(0,nbins)]

    dict_main = {'snps': [num_snps], 'snps_in_annot': [num_sel_snps], 
                 'N': [num_indivs], 'N_covariates': [num_covariates]}
    if detailed:
        dict_detail = {'bin_' + str(idx) + '_N_SNPs': [snps_per_bin[idx]] for idx in range(0, nbins)}
        dict_detail.update({'N_SNPs_per_bin': [num_snps_per_block]})
        dict_main.update(dict_detail)

    return nbins, pd.DataFrame(dict_main).reset_index(drop=True)


def make_base_table(pheno, ancestry):
    return pd.DataFrame({'phenotype_id': [pheno], 'ancestry': [ancestry]}).reset_index(drop=True)


def filter_to_output(list_of_lines):
    match_id = _get_first_idx(list_of_lines, r'OUTPUT:')
    return list_of_lines[match_id:(len(list_of_lines))]


def obtain_variance_and_bins(list_of_lines):
    match_id = _get_first_idx(list_of_lines, r'Variances:')
    post_var_lines = list_of_lines[(match_id+1):(len(list_of_lines))]
    bins = [f for f in post_var_lines if re.search(r'^Sigma\^2_[0-9]{1,}:', f)]
    nbins = len(bins)  
    var_estimate = [_get_first_match(post_var_lines, r'^Sigma\^2_' + str(idx) + r': (' + SCINOT + r'{1,})', transform=float) 
                    for idx in range(0,nbins)]
    var_e = _get_first_match(post_var_lines, r'^Sigma\^2_e: (' + SCINOT + r'{1,})', transform=float)

    dct_var = {'variance_bin_' + str(idx): [var_estimate[idx]] for idx in range(0,nbins)}
    dct_var.update({'variance_e': [var_e]})

    return pd.DataFrame(dct_var).reset_index(drop=True), nbins


def obtain_overlapping_h2(list_of_lines, nbins:int):
    match_id = _get_first_idx(list_of_lines, r'^h\^2\'s \(heritabilities\).+Equation 9 \(overlapping.+$')
    list_of_lines_trim = list_of_lines[match_id:len(list_of_lines)-1]

    # get h2
    regex = r'^.+bin [0-9]{1,} : ('+ SCINOT + r'{1,}) ,  se: (' + SCINOT + r'{1,})'
    match_id = _get_first_idx(list_of_lines_trim, r'^h\^2\'s:')
    list_of_lines_trim = list_of_lines_trim[match_id:len(list_of_lines_trim)]
    vals = _get_two_matches_per_lines(list_of_lines_trim[1:nbins+1], regex, float)
    h2_se = {'h2_s_overlapping_bin_' + str(idx): [val] for idx, (val,_) in enumerate(vals)}
    h2_se.update({'h2_s_se_overlapping_bin_' + str(idx): [se] for idx, (_,se) in enumerate(vals)})

    # get frac
    regex = r'^.+bin [0-9]{1,} : ('+ SCINOT + r'{1,}),  se: ('+ SCINOT + r'{1,})'
    match_id = _get_first_idx(list_of_lines_trim, r'^h\^2_i/h\^2_t:')
    list_of_lines_trim = list_of_lines_trim[match_id:len(list_of_lines_trim)]
    vals = _get_two_matches_per_lines(list_of_lines_trim[1:nbins+1], regex, float)
    h2_frac_se = {'h2_i_h2_t_overlapping_bin_' + str(idx): [val] for idx, (val,_) in enumerate(vals)}
    h2_frac_se.update({'h2_i_h2_t_se_overlapping_bin_' + str(idx): [se] for idx, (_,se) in enumerate(vals)})

    # get enrichments
    regex = r'^.+bin [0-9]{1,} : ('+ SCINOT + r'{1,}),  se: ('+ SCINOT + r'{1,})'
    match_id = _get_first_idx(list_of_lines_trim, r'^Enrichments:')
    list_of_lines_trim = list_of_lines_trim[match_id:len(list_of_lines_trim)]
    vals = _get_two_matches_per_lines(list_of_lines_trim[1:nbins+1], regex, float)
    enr_se = {'Enrich_overlapping_bin_' + str(idx): [val] for idx, (val,_) in enumerate(vals)}
    enr_se.update({'Enrich_se_overlapping_bin_' + str(idx): [se] for idx, (_,se) in enumerate(vals)})

    final_dict = h2_se
    final_dict.update(h2_frac_se)
    final_dict.update(enr_se)

    return pd.DataFrame(final_dict).reset_index(drop=True)


def obtain_nonoverlapping_h2(list_of_lines, nbins:int):
    match_id = _get_first_idx(list_of_lines, r'^h\^2\'s \(heritabilities\).+Equations 2\-4 \(non\-overlapping.+$')
    list_of_lines_trim = list_of_lines[match_id:len(list_of_lines)]

    # get h2
    regex = r'^.+bin [0-9]{1,} : ('+ SCINOT + r'{1,}),  se: ('+ SCINOT + r'{1,})'
    match_id = _get_first_idx(list_of_lines_trim, r'^h\^2\'s:')
    list_of_lines_trim = list_of_lines_trim[match_id:len(list_of_lines_trim)]
    vals = _get_two_matches_per_lines(list_of_lines_trim[1:nbins+1], regex, float)
    h2_se = {'h2_s_non_overlapping_bin_' + str(idx): [val] for idx, (val,_) in enumerate(vals)}
    h2_se.update({'h2_s_se_non_overlapping_bin_' + str(idx): [se] for idx, (_,se) in enumerate(vals)})

    # get h2_total
    regex = r'^Total h\^2 : (' + SCINOT + r'{1,}), se: (' + SCINOT + r'{1,})'
    match_id = _get_first_idx(list_of_lines_trim, regex)
    list_of_lines_trim = list_of_lines_trim[match_id:len(list_of_lines_trim)]
    vals = _get_two_matches_per_lines([list_of_lines_trim[0]], regex, float)[0]
    h2_tot = {'h2': vals[0], 'h2_se': vals[1]}

    # get enrichments
    regex = r'^.+bin [0-9]{1,}: (' + SCINOT + r'{1,}) ,  se: (' + SCINOT + r'{1,})'
    match_id = _get_first_idx(list_of_lines_trim, r'^Enrichments:')
    list_of_lines_trim = list_of_lines_trim[match_id:len(list_of_lines_trim)]
    vals = _get_two_matches_per_lines(list_of_lines_trim[1:nbins+1], regex, float)
    enr_se = {'Enrich_non_overlapping_bin_' + str(idx): [val] for idx, (val,_) in enumerate(vals)}
    enr_se.update({'Enrich_se_non_overlapping_bin_' + str(idx): [se] for idx, (_,se) in enumerate(vals)})

    final_dict = h2_tot
    final_dict.update(h2_se)
    final_dict.update(enr_se)

    return pd.DataFrame(final_dict).reset_index(drop=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', type=str,
                        help='RHE-mc (HE) log file to parse. This script will intuit the number of ' + \
                            'bins analyzed.')
    parser.add_argument('--pheno', type=str,
                        help='Full phenotype code analyzed for reference purposes.')
    parser.add_argument('--ancestry', type=str,
                        help='Ancestry group analyzed for reference purposes.')
    parser.add_argument('--include-trivial-details', action='store_true',
                        help='Include very specific details including the number of SNPs in each bin and overlapping heritability estimates.')
    parser.add_argument('--out', type=str,
                        help='Output location for parsed table. Expects this to contain the extension. Output ' + \
                            'is tab-delimited with a single row, but with included column headers.')

    args = parser.parse_args()

    list_of_lines = filter_to_output(read_file(args.file))
    this_table = make_base_table(args.pheno, args.ancestry)

    variance, nbins = obtain_variance_and_bins(list_of_lines)

    overlapping = obtain_overlapping_h2(list_of_lines, nbins)
    non_overlapping = obtain_nonoverlapping_h2(list_of_lines, nbins)

    if args.include_trivial_details:
        tab = pd.concat([this_table, non_overlapping, overlapping, variance], axis=1)
    else:
        tab = pd.concat([this_table, non_overlapping, variance], axis=1)

    tab.to_csv(path_or_buf=args.out, header=True, index=False, sep='\t', na_rep='NA')