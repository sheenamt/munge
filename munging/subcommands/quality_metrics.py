"""
Parse picard rmdup and hs metrics output to create quality metrics file

Usage:

munge quality_metrics $SAVEPATH/$PFX.quality_metrics $SAVEPATH/$PFX.hs_metrics -o $SAVEPATH/$PFX_quality_metrics_Analysis.txt
"""
import argparse
import csv
import sys
import re
from numpy import array, average

import pprint

def build_parser(parser):
    parser.add_argument(
        'metric_file', type=argparse.FileType('rU'),
        default=sys.stdin)
    parser.add_argument(
        'hsmetrics', type=argparse.FileType('rU'),
        default=sys.stdin)
    parser.add_argument(
        '-o', '--outfile',
        help='Output file', default=sys.stdout,
        type=argparse.FileType('w'))

def action(args):

    metric_info = csv.reader(args.metric_file, delimiter="\t")
    hs_info = csv.reader(args.hsmetrics, delimiter='\t')
    writer = csv.writer(args.outfile, quoting=csv.QUOTE_MINIMAL, delimiter='\t')
    for line in metric_info:
        try:
            if re.search('LIBRARY', line[0]):
                writer.writerow(line)
                writer.writerow(next(metric_info))
        except IndexError:
            continue

    #TOTAL_READS:5, PF_UNIQUE_READS:7, PF_UQ_READS_ALIGNED:10, PCT_SELECTED_BASES:17, PCT_OFF_BAIT:18,
    #MEAN_TARGET_COVERAGE:21, PCT_USABLE_BASES_ON_TARGET:23, ZERO_CVG_TARGETS_PCT:25, AT_DROPOUT:35, GC_DROPOUT:36, LIBRARY:38
    items = [5, 7, 10, 17, 18, 21, 23, 25, 35, 36]
    new_line = []
    for line in hs_info:
        try:
            if not re.search('^#', line[0]):
                for i, value in enumerate(line):
                    if i in items:
                        new_line.append(value)
                writer.writerow(new_line)
                #reset new_line so the header isn't printed twice
                new_line = []

        except IndexError:
            continue

