#!/usr/bin/env python

__author__ = 'Rahul Gupta'

import re
import argparse
from parse_ldsc_log import import_log

def internal_parse_args(parser):
    """ Parse incoming arguments, catch issues, and convert to lists.

    Keyword arguments:

    parser -- an ArgumentParser()

    Returns:
    dict of arguments: new values
    """
    args = parser.parse_args()
    if args.log is None:
        raise ValueError('Must provide at one log file.')
    return args


parser = argparse.ArgumentParser()
parser.add_argument('--log', type=str, default=None,
                    help='Single log to test.')

if __name__ == '__main__':
    args = internal_parse_args(parser)
    log = import_log(args.log)
    stop_str = '5 \(or 6 in the case of case\/control\) ancestry specific columns not found; exiting.'
    tf_skip = False
    
    for l in log:
        if re.search(stop_str, l):
            tf_skip = True
            break
        else:
            continue
    
    print('Incomplete' if tf_skip else 'Complete')