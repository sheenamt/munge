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
import pandas as pd
from munging import filters
from munging.utils import walker, munge_pfx
from natsort import natsorted

def build_parser(parser):
    parser.add_argument('path',
                        help='Path to analysis files')
    parser.add_argument('pipeline_manifest', type=argparse.FileType('rU'),
                        help='Path to pipeline manifest, used for ordering output')    
    parser.add_argument('outfile', type = argparse.FileType('w'),
                        default = sys.stdout,
                        help='Name of the output file')
    

def action(args):
    #Grab all analysis files from the path 
    files = ifilter(filters.any_analysis, walker(args.path))        
    files = filter(filters.polyhunter_analysis, files)
    #    interesting_files = glob.glob("*.csv")
    df_list = []
    df = pd.DataFrame()
    sort_order = [x['barcode_id'] for x in csv.DictReader(args.pipeline_manifest)]
    for sample in sort_order:
        #Grab the file for each sample, in specified sort order
        pfx_file = [s for s in files if sample in s.fname]
        if pfx_file:
            pfx_file = pfx_file[0]
            pfx = munge_pfx(pfx_file.fname)
            #Create a smaller version of this really long string
            data = pd.read_csv(os.path.join(pfx_file.dir, pfx_file.fname), sep='\t')
            data.index = [pfx['mini-pfx']] * len(data)
            df = df.append(data)
    cols=natsorted(df.columns)
    df.to_csv(args.outfile, sep='\t',na_rep='0',columns=cols)

