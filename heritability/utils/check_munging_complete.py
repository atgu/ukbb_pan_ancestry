#!/usr/bin/env python

__author__ = 'Rahul Gupta'

import re
import argparse
from parse_ldsc_log import import_log


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--log', type=str, default=None, required=True,
                        help='Single log to test.')

    args = parser.parse_args()
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