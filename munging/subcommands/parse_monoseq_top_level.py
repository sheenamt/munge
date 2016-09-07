"""Parse the polyhunter output, give info for all variants
"""
import os
import csv
import re
import sys
import copy
import argparse
import collections
from itertools import count, groupby, chain, ifilter , izip_longest
from operator import itemgetter

from munging import filters
from munging.utils import walker, munge_pfx

def build_parser(parser):
    parser.add_argument('path',
                        help='Path to analysis files')
    parser.add_argument('-o','--outfile', type = argparse.FileType('w'),
                        default = sys.stdout,
                        help='Name of the output file')
    

def action(args):
    specimens = collections.defaultdict(dict)
    #Grab all analysis files from the path 
    files = ifilter(filters.any_analysis, walker(args.path))        
    files = ifilter(filters.polyhunter_analysis, files)
    #sort the files so that the output in the workbook is sorted 
    files=sorted(files)
    annotation_header = set(['Sample',])
    for pth in files:
        pfx = munge_pfx(pth.fname)
        with open(os.path.join(pth.dir, pth.fname)) as fname:
            reader = csv.DictReader(fname, delimiter='\t')            
            for key in reader.fieldnames:
                annotation_header.add(key)
            for row in reader:
                specimens[pfx['mini-pfx']]=row

    fieldnames = sorted(annotation_header, reverse=True)
    writer = csv.DictWriter(args.outfile, fieldnames = fieldnames,  extrasaction = 'ignore', delimiter = '\t')
    writer.writeheader()
    for key, value in specimens.items():
        row = {'Sample':key}
        row.update(value)
        writer.writerow(row)
