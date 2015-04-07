"""
Create the summary file of your choice

Usage:

munge create_summary <Type> $SAVEPATH -o $OUTFILE

"""
import subprocess
import tempfile
import logging
import shutil
import os
import re
import csv
import sys
import argparse
import sqlite3
from collections import defaultdict, namedtuple
import collections
from itertools import count, groupby, chain, ifilter

from munging.utils import dict_factory
from munging.annotation import split_chr_loc, multi_split
from munging import filters
from munging import parsers 
from munging.utils import walker, munge_pfx

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


