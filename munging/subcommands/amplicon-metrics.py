"""
Parse Illumina Amplicon output to create amplicon quality file

Usage:

munge data/HHv1.csv amplicon_metrics_file $SAVEPATH/ -o $SAVEPATH/$PFX_amplicon_metrics_Analysis.txt
"""
import argparse
import sys
import csv

from munging import parsers

def build_parser(parser):
    parser.add_argument(
        'amplicons', type=argparse.FileType('rU'),
        default=sys.stdin)
    parser.add_argument(
        'run_metrics', type=argparse.FileType('rU'),
        default=sys.stdin)
    parser.add_argument(
        '-o', '--outfile',
        help='Output file', default=sys.stdout,
        type=argparse.FileType('w'))

def action(args):

    run_metrics=args.run_metrics
    amplicons=args.amplicons
    import IPython
    IPython.embed()
    lines=metricsfile.readlines()
    #filter out lines that start with # or are just new lines
    #then a simplelines of code to make a dictionary and then print
    variant_keys =[]
    output_dict, variant_keys = parsers.parse_hsmetrics(lines, variant_keys)
    writer = csv.DictWriter(args.outfile, fieldnames=variant_keys, quoting=csv.QUOTE_MINIMAL,extrasaction='ignore', delimiter='\t')
    writer.writeheader()
    writer.writerow(output_dict)
