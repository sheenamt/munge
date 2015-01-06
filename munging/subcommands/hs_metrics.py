"""
Parse picard and CNV output to create quality metrics file

Usage:

munge hs_metrics $SAVEPATH/$PFX.hs_metrics -o $SAVEPATH/$PFX_quality_metrics_Analysis.txt
"""
#!/usr/bin/python
import argparse
import sys
import os
import string
import re
import IPython
from munging.annotation import multi_split
import itertools
import csv
""" parse a picard HsMetrics *.metrics output to something more human readable """

def build_parser(parser):
    parser.add_argument(
        'metricsfile', type=argparse.FileType('rU'),
        default=sys.stdin)
    parser.add_argument(
        '-o', '--outfile',
        help='Output file', default=sys.stdout,
        type=argparse.FileType('w'))

def action(args):
    metricsfile=args.metricsfile
    samplename=multi_split(metricsfile.name,'/.')[-2]
    lines=metricsfile.readlines()
    #filter out lines that start with # or are just new lines
    datalines=filter(lambda x: '#' not in x, lines)
    datalines=filter(lambda x: x!='\n' , datalines)
    #then strip newlines from the lines that contain data
    datalines=map(lambda x: x.strip() , datalines)
    #then a simplelines of code to make a dictionary and then print
    keys=datalines[0].split('\t')
    values=datalines[1].split('\t')
    #metrics_dict=dict(zip(keys,values))
    metrics_dict = dict(itertools.izip_longest(keys,values, fillvalue='NA'))
    writer = csv.writer(args.outfile, quoting=csv.QUOTE_MINIMAL, delimiter='\t')
    #Only print the keys we care about:
    print_keys=['TOTAL_READS',
                'PF_UNIQUE_READS',
                'PF_UQ_READS_ALIGNED',
                'PCT_SELECTED_BASES',
                'PCT_OFF_BAIT',
                'MEAN_TARGET_COVERAGE',
                'PCT_USABLE_BASES_ON_TARGET',
                'ZERO_CVG_TARGETS_PCT',
                'AT_DROPOUT',
                'GC_DROPOUT',
                'LIBRARY']

    for k in print_keys:
        writer.writerow([ k, metrics_dict[k], samplename ])

