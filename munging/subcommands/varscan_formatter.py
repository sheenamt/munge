"""
Parses varscan readcounts for specific positions 
"""
import sys
import subprocess
import tempfile
import logging
import shutil
import os
import pandas as pd

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('clin_flagged', 
        help='A required input file')
    parser.add_argument('outfile', 
        help='A required output file')

def fix_deletions(series):
    if '-' in series['var_base']:
        return pd.Series([series['chr'],
                          series['start']-1,
                          '-',
                          '-'.join(['DEL',str(len(series['ref_base'])),series['ref_base']])],
                         index=['chr','start','ref_base','var_base'])
    elif '-' in series['ref_base']:
        return pd.Series([series['chr'],
                          series['start'],
                          '-',
                          '-'.join(['INS',str(len(series['var_base'])),series['var_base']])],
                         index=['chr','start','ref_base','var_base'])
    else:
        return pd.Series([series['chr'],series['start'],series['ref_base'],series['var_base']],
                         index=['chr','start','ref_base','var_base'])

def action(args):

    flagged_variants=args.clin_flagged

    cols=['chr','start','stop','ref_base','var_base','comment']
    variants=pd.read_csv(flagged_variants, delimiter='\t', header=None, index_col=False, names = cols)
    genotype_calls=args.outfile
    
    variants.apply(fix_deletions, axis = 1).to_csv(genotype_calls, index=False, header=None,sep='\t')
    
    
