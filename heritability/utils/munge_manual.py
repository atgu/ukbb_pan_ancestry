#!/usr/bin/env python

__author__ = 'Rahul Gupta'

import re, logging, argparse
import pandas as pd
import numpy as np
from scipy.stats import chi2
from itertools import compress


def read_file(filename):
    """ Read summary statistis file using Pandas. Exports import outcome to log.

    Parameters
    ----------
    filename : :obj: `str`

    Returns
    -------
    :obj: `DataFrame`
        Summary statistics.
    """
    data_tsv = pd.read_csv(filename,sep="\t")
    logging.info('File %s imported successfully.', filename)
    return data_tsv


def search_column(table, string, enforce_singular = True):
    """ Searches column names for strings.

    Parameters
    ----------
    table : :obj: `DataFrame`

    string : :obj: `str`
        String to search column names for.

    enforce_singular : :obj: `bool`
        If true, will throw a ValueError if there are 0 or more than 1 matches
        with the provided string. 

    Returns
    -------
    :obj: `DataFrame`
        Summary statistics.
    """
    column_names = list(table.columns)
    found_list = list(compress(column_names, [re.search(string,col) for col in column_names]))
    if (len(found_list) != 1) & (enforce_singular): 
        raise ValueError('There must exactly be one ' + string + ' column in the inputted summary statistics. Currently found ' + str(len(found_list)) + '.')
    return(found_list[0])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--sumstats', type=str,
                        help="Filename to munge.")
    parser.add_argument('--N', type=float,
                        help="Sample size. We assume that each SNP has the same N.")
    parser.add_argument('--out', type=str,
                        help="Output filename. Will output as tab-delimited.")
    parser.add_argument('--logfile', type=str,
                        help="Name of log file to output.")
    
    args = parser.parse_args()
    if args.sumstats is None:
        raise ValueError('--sumstats must be provided as a path to the extracted sumstats file.')
    if args.N is None:
        raise ValueError('--N must be provided.')
    if args.out is None:
        raise ValueError('--out must be provided. This is the output location and filename.')
    if args.logfile is None:
        raise ValueError('--logifle must be provided. This is the output location of the logfile.')

    logging.basicConfig(filename = args.logfile, level = logging.DEBUG)

    data = read_file(args.sumstats)
    data["SNP"] = data["chr"].astype('str') + ":" + data["pos"].astype('str') + ":" + data["ref"] + ":" + data["alt"]
    data["A1"] = data["ref"]
    data["A2"] = data["alt"]
    data["N"] = args.N
    logging.info('Used N = %s.', str(args.N))
    logging.info('Renamed ref -> A1 and alt -> A2.')
    logging.info('Concatenated chr, pos, ref, alt -> SNP')

    pval_col_touse = search_column(data, 'pval', enforce_singular = True)
    logging.info('Using %s for pval.', pval_col_touse)
    beta_col_touse = search_column(data, 'beta', enforce_singular = True)
    logging.info('Using %s for beta.', beta_col_touse)

    # the af field is either af or af_controls/cases for continuous and
    # categorical variables respectively.
    try:
        af_col_touse = search_column(data, 'af', enforce_singular = True)
    except ValueError as err:
        logging.info('The following arose for af: ' + str(err[0]))
        logging.info('Trying af_controls.')
        af_col_touse = search_column(data, 'af_controls', enforce_singular = True)
    logging.info('Using %s for af.', af_col_touse)

    low_conf_col_touse = search_column(data, 'low_confidence', enforce_singular = True)
    logging.info('Using %s for low_confidence.', low_conf_col_touse)
    logging.info('-----------------')

    original_row_count = data.shape[0]
    logging.info('Original row count: %d', original_row_count)
    data_nona = data.loc[~(data[pval_col_touse].isnull() | \
                           data[beta_col_touse].isnull() | \
                           data[af_col_touse].isnull() | \
                           data[low_conf_col_touse].isnull())]
    row_count_after_na_remove = data_nona.shape[0]
    logging.info('Rows removed due to NA in %s, %s, %s, or %s: %d',
                 pval_col_touse,
                 beta_col_touse, 
                 af_col_touse,
                 low_conf_col_touse,
                 original_row_count - row_count_after_na_remove)

    # filter to p value column within bounds
    data_p_filt = data_nona.loc[(data_nona[pval_col_touse] > 0) & \
                                (data_nona[pval_col_touse] <= 1)]
    post_p_filt_nrow = data_p_filt.shape[0]
    logging.info('Rows removed due to %s out of bounds: %d',
                 pval_col_touse,
                 row_count_after_na_remove - post_p_filt_nrow)

    # filter to af >= 0.01
    data_af_filt = data_p_filt.loc[(data_p_filt[af_col_touse] >= 0.01) & \
                                   (data_p_filt[af_col_touse] <= 1)]
    post_af_filt_nrow = data_af_filt.shape[0]
    logging.info('Rows removed due to %s below 0.01: %d',
                 af_col_touse,
                 post_p_filt_nrow - post_af_filt_nrow)

    # remove low confidence vars
    data_conf_filt = data_af_filt.loc[data_af_filt[low_conf_col_touse].astype("bool") == False]
    post_lowconf_filt_nrow = data_conf_filt.shape[0]
    logging.info('Rows removed due to low confidence: %d',
                 post_af_filt_nrow - post_lowconf_filt_nrow)

    # compute z
    data_conf_filt.loc[:, "Z"] = np.sqrt(chi2.isf(data_conf_filt[pval_col_touse], 1))

    # attach sign to z based on the sign of beta
    data_conf_filt.loc[:, "Z"] *= np.sign(data_conf_filt[beta_col_touse])

    # obtain final result
    data_final = data_conf_filt.loc[:, [
        "SNP", "A1", "A2", "N", "Z", beta_col_touse, pval_col_touse, af_col_touse]]

    # provide metrics as in munge_sumstats from ldsc
    logging.info('Final number of rows in file: %d', data_final.shape[0])
    logging.info('-----------------')
    logging.info('\nMetadata:')
    chisquare_values = (data_final.Z ** 2)
    mean_chisq = chisquare_values.mean()
    logging.info('Mean chi^2 = ' + str(round(mean_chisq, 3)))
    if mean_chisq < 1.02: 
        logging.warning("WARNING: mean chi^2 may be too small.")
    logging.info('Lambda GC = ' + str(round(chisquare_values.median() / 0.4549, 3)))
    logging.info('Max chi^2 = ' + str(round(chisquare_values.max(), 3)))
    logging.info('%d genome-wide significant SNPs (some may have been removed by filtering)', 
                 (chisquare_values > 29).sum())

    # output
    data_final.to_csv(args.out, sep='\t', compression='gzip', index=False, float_format='%.3f')
