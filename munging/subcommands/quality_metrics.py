"""
Parse picard rmdup and hs metrics output to create quality metrics file

Usage:

munge quality_metrics $SAVEPATH/$PFX.quality_metrics $SAVEPATH/$PFX.hs_metrics -o $SAVEPATH/$PFX_quality_metrics_Analysis.txt
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
        'hsmetrics', type=argparse.FileType('rU'),
        default=sys.stdin)
    parser.add_argument(
        '-o', '--outfile',
        help='Output file', default=sys.stdout,
        type=argparse.FileType('w'))

def action(args):
    variant_keys = []
    quality_metricsfile=args.quality_metrics
    lines=quality_metricsfile.readlines()
    quality_dict, quality_keys = parsers.parse_qualitymetrics(lines, variant_keys)    

    metricsfile=args.hsmetrics
    lines=metricsfile.readlines()
    hsmetrics_dict, hsmetrics_keys = parsers.parse_hsmetrics(lines, variant_keys)

    variant_keys = quality_keys + hsmetrics_keys
    output_dict = dict(quality_dict,**hsmetrics_dict)
    writer = csv.DictWriter(args.outfile, fieldnames=variant_keys, quoting=csv.QUOTE_MINIMAL,extrasaction='ignore', delimiter='\t')
    writer.writeheader()
    writer.writerow(output_dict)
