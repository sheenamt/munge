"""
Parse picard hs-metrics output to create quality metrics file

Usage:

munge hs_metrics $SAVEPATH/$PFX.hs_metrics -o $SAVEPATH/$PFX_quality_metrics_Analysis.txt
"""
import argparse
import sys
import csv

from munging import parsers

def build_parser(parser):
    parser.add_argument(
        'hsmetrics', type=argparse.FileType('rU'),
        default=sys.stdin)
    parser.add_argument(
        '-o', '--outfile',
        help='Output file', default=sys.stdout,
        type=argparse.FileType('w'))

def action(args):
    metricsfile=args.hsmetrics
    lines=metricsfile.readlines()
    #filter out lines that start with # or are just new lines
    #then a simplelines of code to make a dictionary and then print
    variant_keys =[]
    output_dict, variant_keys = parsers.parse_hsmetrics(lines, variant_keys)
    writer = csv.DictWriter(args.outfile, fieldnames=variant_keys, quoting=csv.QUOTE_MINIMAL,extrasaction='ignore', delimiter='\t')
    writer.writeheader()
    writer.writerow(output_dict)
