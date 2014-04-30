"""
Crawl analysis files to create one analysis file with all info

Usage:

 munge combined_cnv $SAVEPATH -o $SAVEPATH/${RUN}_Combined_CNV

"""

import logging
import os
import csv
import sys
import argparse
import collections
from itertools import count, groupby, chain, ifilter

from munging.utils import dict_factory
from munging import filters
from munging.utils import walker


log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('type', choices = ['Exon','Gene'],
                        help = 'Type of cnv file to parse')
    parser.add_argument('path',
                        help='Path to analysis files')
    parser.add_argument('-o','--outfile', type = argparse.FileType('w'),
                        default = sys.stdout,
                        help='Name of the output file')


def action(args):
    specimens= collections.defaultdict(dict)
    annotation = {}
    prefixes = []
    # apply a series of filters to files
    files = ifilter(filters.any_analysis, walker(args.path))
    if args.type =='Exon':
        files = ifilter(filters.cnv_exon_analysis, files)
    elif args.type=='Gene':
        files = ifilter(filters.cnv_gene_analysis, files)
    variant_keys = ['Position', 'Gene' ]

    #sort the files so that the output in the workbook is sorted
    files=sorted(files)
    for pth in files:
        print pth
        pfx = pth.fname.split('_')[0]
        log.warning('%s %s %s' % (pfx, pth.dir, pth.fname))
        log_pfx=pfx+'_Log'
        prefixes.append(log_pfx)
        with open(os.path.join(args.path, pth.fname)) as fname:
            reader = csv.DictReader(fname, delimiter='\t')
            for row in reader:
                variant = tuple(row[k] for k in variant_keys)
                specimens[variant][log_pfx] = row['Ave_Adjusted_Log_Ratio']
                annotation[variant] = row


    annotation_headers = [
        'Transcripts']

    writer = csv.DictWriter(args.outfile, fieldnames = variant_keys + annotation_headers + prefixes,  extrasaction = 'ignore', delimiter='\t')
    writer.writeheader()
    for variant in sorted(specimens.keys()):
        d = {k:v for k,v in zip(variant_keys,variant)}
        d.update({pfx:specimens[variant].get(pfx) for pfx in prefixes})
        d.update(annotation[variant])
        writer.writerow(d)
