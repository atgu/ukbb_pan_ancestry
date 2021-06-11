#!/usr/bin/env python

__author__ = 'Rahul Gupta'

import argparse
import pandas as pd


def internal_parse_args(parser):
    """ Parse incoming arguments, catch issues, and convert to lists.

    Keyword arguments:

    parser -- an ArgumentParser()

    Returns:
    dict of arguments: new values
    """
    args = parser.parse_args()
    if args.tables is None:
        raise ValueError('Must provide at least one table file.')
    if args.out is None:
        raise ValueError('Must provide output directory.')
    if args.read_list_disk:
        list_of_tables = list(pd.read_csv(args.tables, sep='\t', names=['gs']).gs)
    else:
        list_of_tables = args.tables.split(',')
    args_dict = vars(args)
    args_dict.update({'tables': list_of_tables})
    return args_dict


parser = argparse.ArgumentParser()
parser.add_argument('--tables', type=str, default=None,
                    help='Comma delimited list of tables to import.')
parser.add_argument('--read-list-disk', action='store_true',
                    help='If provided, will assume that the tables argument points ' + \
                    'to a local file in which each line is a path to a table.')
parser.add_argument('--out', type=str, default=None,
                    help='Output directory.')

if __name__ == '__main__':
    args = internal_parse_args(parser)
    tabs = [pd.read_csv(tab, sep='\t') for tab in args['tables']]
    final_table = pd.concat(tabs)
    final_table.to_csv(sep='\t', path_or_buf=args['out'], index=False)