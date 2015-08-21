"""
Create the summary file of your choice

Usage:

munge create_summary <Type> $SAVEPATH -o $OUTFILE

"""
import logging
import csv
import sys
import argparse
import collections
from itertools import ifilter

from munging import filters, parsers
from munging.utils import walker

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('type', 
                        choices=['pindel','snp','cnv_exon','cnv_gene','quality','msi','clin_flagged','msi_flagged'],
                        help='Type of output summary to create')
    parser.add_argument('path',
                        help='Path to analysis files')
    parser.add_argument('-o','--outfile', type = argparse.FileType('w'),
                        default = sys.stdout,
                        help='Name of the output file')
    

def action(args):
    specimens = collections.defaultdict(dict)
    annotation = {}
    prefixes = []
    variant_keys = []
    #Grab all analysis files from the path 
    files = ifilter(filters.any_analysis, walker(args.path))        
    analysis_type='_'.join(['parsers.parse',args.type])
    print "analysis type:",analysis_type
    chosen_parser='{}(files, specimens, annotation, prefixes, variant_keys)'.format(analysis_type)
    specimens, annotation, prefixes, fieldnames, variant_keys=eval(chosen_parser)
        

    writer = csv.DictWriter(args.outfile, fieldnames = fieldnames,  extrasaction = 'ignore', delimiter = '\t')
    writer.writeheader()
    for variant in sorted(specimens.keys()):
        d = {k:v for k,v in zip(variant_keys,variant)}
        d.update({pfx:specimens[variant].get(pfx) for pfx in prefixes})
        d.update(annotation[variant])
        writer.writerow(d)


