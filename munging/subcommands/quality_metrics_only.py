"""
Parse picard rmdup output to create quality metrics file

Usage:

munge quality_metrics_only $SAVEPATH/$PFX.quality_metrics -o $SAVEPATH/$PFX_quality_metrics_Analysis.txt
"""
import argparse
import csv
import sys

from munging import parsers

def build_parser(parser):
    parser.add_argument(
        'quality_metrics', type=argparse.FileType('rU'),
        default=sys.stdin)
    parser.add_argument(
        '-o', '--outfile',
        help='Output file', default=sys.stdout,
        type=argparse.FileType('w'))

def action(args):

    quality_metricsfile=args.quality_metrics
    lines=quality_metricsfile.readlines()
    #filter out lines that start with # or are just new lines
    #then a simplelines of code to make a dictionary and then print
    variant_keys=[]
    output_dict, variant_keys = parsers.parse_qualitymetrics(lines, variant_keys)    
    writer = csv.DictWriter(args.outfile, fieldnames=variant_keys, quoting=csv.QUOTE_MINIMAL,extrasaction='ignore', delimiter='\t')
    writer.writeheader()
    writer.writerow(output_dict)
