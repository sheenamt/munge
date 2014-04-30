"""
Crawl analysis files to create one analysis file with all info

Usage:

 munge combined_pindel $SAVEPATH -o $SAVEPATH/${RUN}_Combined_Pindel

"""

import logging
import os
import csv 
import sys
import argparse 
import collections
from itertools import count, groupby, chain, ifilter

from munging.utils import dict_factory
from munging.annotation import split_chr_loc
from munging import filters
from munging.utils import walker

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('path', 
                        help='Path to analysis files')
    parser.add_argument('-o','--outfile', type = argparse.FileType('w'),
                        default = sys.stdout,
                        help='Name of the output file')


def parse_pindel(variant_keys, files, path):
    specimens = collections.defaultdict(dict)
    annotation = {}
    prefixes = []

    for pth in files:
        pfx = pth.fname.split('_')[0]
        log.warning('%s %s %s' % (pfx, pth.dir, pth.fname))
        prefixes.append(pfx)
        with open(os.path.join(path, pth.fname)) as fname:
            print pth.fname
            reader = csv.DictReader(fname, delimiter='\t')
            for row in reader:
                variant = tuple(row[k] for k in variant_keys)
                specimens[variant][pfx] = row['Reads']
                annotation[variant] = row
    return specimens, annotation, prefixes

def action(args):

    files = ifilter(filters.any_analysis, walker(args.path))
    files = ifilter(filters.pindel_analysis, files)
    #sort the files so that the output in the workbook is sorted
    files=sorted(files)

    variant_keys = ['Position', 'Gene']
    specimens, annotation, prefixes=parse_pindel(variant_keys, files, args.path)
    annotation_headers = [
        'Gene_Region',
        'Event_Type',
        'Size',
        'Transcripts',
        ]
                        
    writer = csv.DictWriter(args.outfile, fieldnames = variant_keys + annotation_headers + prefixes,  extrasaction = 'ignore', delimiter = '\t')    
    writer.writeheader()
    for variant in sorted(specimens.keys()):                
        d = {k:v for k,v in zip(variant_keys,variant)}
        d.update({pfx:specimens[variant].get(pfx) for pfx in prefixes})
        d.update(annotation[variant])
        writer.writerow(d)
