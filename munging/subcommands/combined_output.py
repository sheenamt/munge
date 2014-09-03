"""
Crawl analysis files to create one analysis file with all info

Usage:

 munge combined_output $SAVEPATH -o $SAVEPATH/${RUN}_Combined_Analysis


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
import pprint

from munging.utils import dict_factory
from munging.annotation import split_chr_loc, multi_split
from munging import filters
from munging.utils import walker

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument(
        'infiles', action='append', nargs='+',
        help='Input files')
    parser.add_argument('-o','--outfile', type = argparse.FileType('w'),
                        default = sys.stdout,
                        help='Name of the output file')


def action(args):
    specimens = collections.defaultdict(dict)
    annotation = {}
    prefixes = []
    variant_keys = ['Position', 'Ref_Base', 'Var_Base']
    (infiles, ) = args.infiles
    print "infiles:", infiles
    files = ifilter(filters.any_analysis, infiles)
    print "analysis files:", files
    files = ifilter(filters.only_analysis, files)
    #sort the files so that the output in the workbook is sorted
    files=sorted(files)

    for pth in files:
        print pth
        pfx = pth.split('_')[0]
        # ref_pfx=pfx+'_Ref'
        # var_pfx=pfx+'_Var'
        reads_pfx=pfx+'_Ref|Var'
        # prefixes.append(ref_pfx)
        # prefixes.append(var_pfx)
        prefixes.append(reads_pfx)
        with open(os.path.join(pth)) as fname:
            reader = csv.DictReader(fname, delimiter='\t')
            for row in reader:
                variant = tuple(row[k] for k in variant_keys)
                # specimens[variant][ref_pfx] = row['Ref_Reads']
                # specimens[variant][var_pfx] = row['Var_Reads']
                specimens[variant][reads_pfx] = row['Ref_Reads']+'|'+row['Var_Reads']
                annotation[variant] = row


    annotation_headers = [
        'Gene',
        'Variant_Type',
        'Transcripts',
        'Clinically_Flagged',
        'Cosmic',
        'Segdup',
        'Polyphen',
        'Sift',
        'Mutation_Taster',
        'Gerp',
        'HiSeq_Freq',
        'HiSeq_Count',
        'MiSeq_Freq',
        'MiSeq_Count',
        '1000g_ALL',
        'EVS_esp6500_ALL',
        '1000g_AMR',
        'EVS_esp6500_AA',
        '1000g_EUR',
        'EVS_esp6500_EU',
        '1000g_ASN',
        '1000g_AFR']

    writer = csv.DictWriter(args.outfile, fieldnames = variant_keys + annotation_headers + prefixes,  extrasaction = 'ignore', delimiter = '\t')    
    writer.writeheader()
    for variant in sorted(specimens.keys()):                
        d = {k:v for k,v in zip(variant_keys,variant)}
        d.update({pfx:specimens[variant].get(pfx) for pfx in prefixes})
        d.update(annotation[variant])
        writer.writerow(d)

