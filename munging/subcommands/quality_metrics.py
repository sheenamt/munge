"""
Parse picard and CNV output to create quality metrics file

Usage:

munge quality_metrics $SAVEPATH/$PFX.quality_metrics $SAVEPATH/$PFX.hs_metrics -o $SAVEPATH/$PFX_quality_metrics_Analysis.txt
"""
import argparse
import csv
import sys
import re
from numpy import array, average


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


def get_values(info):
    """
    Get the standard deviation from the tumor.rd.ori column in the CNV file
    """
    tumors = []
    for i in info:
        tumors.append(float(i[16]))  # this is the tumor.rd.ori column
    b = array(tumors)
    return b


def action(args):

    metric_info = csv.reader(args.metric_file, delimiter="\t")
    writer = csv.writer(args.outfile, quoting=csv.QUOTE_MINIMAL, delimiter='\t')
    for line in metric_info:
        try:
            if re.search('LIBRARY', line[0]):
                writer.writerow(line)
                writer.writerow(next(metric_info))
        except IndexError:
            continue

    hs_info = csv.reader(args.hsmetrics, delimiter='\t')
    #TOTAL_READS, PF_UNIQUE_READS, PF_UQ_READS_ALIGNED, PCT_SELECTED_BASES, PCT_OFF_BAIT,
    #MEAN_TARGET_COVERAGE, PCT_USABLE_BASES_ON_TARGET, ZERO_CVG_TARGETS_PCT, AT_DROPOUT, GC_DROPOUT, LIBRARY
    items = [5, 7, 10, 17, 18, 21, 23, 25, 35, 36, 38]
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

