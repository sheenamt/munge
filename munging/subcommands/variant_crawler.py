"""
Create annovar file from Clinical variants csv

Usage:

munge variant_crawler /home/genetics/genetics_pipeline/doc/Clinical_variants.csv -o hg19_clinical_variants

"""

import logging
import csv
import argparse
import sys

from munging.annotation import split_chr_loc

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument(
        'infile', type=argparse.FileType('rU'),
        default=sys.stdin)
    parser.add_argument(
        '-o','--outfile',
        help='Output file', default = sys.stdout,
        type = argparse.FileType('w'))

def get_info(fname):
    """
    Function to get information from csv file and create dictionary to write output file
    """
    reader = csv.DictReader(fname, delimiter='\t')
    for d in reader:
         d.update()
         try:
             d['chromosome'], d['start'], d['end']=split_chr_loc(d['chr_loc'])
         except KeyError:
             print 'Key not found', fname
             continue
         yield d



def action(args):
    writer = csv.DictWriter(args.outfile,
                            fieldnames = ['chromosome', 'start', 'end', 'ref_base', 'var_base', 'Clinically_Flagged'],
                            delimiter='\t',
                            extrasaction = 'ignore')
    fname=args.infile
    for d in get_info(fname):
        writer.writerow(d)
