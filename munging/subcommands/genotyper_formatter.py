"""
Parses varscan readcounts for specific positions 
"""
import sys
import subprocess
import logging
import os
import argparse
import pandas as pd
from collections import namedtuple

log = logging.getLogger(__name__)
pd.options.display.width = 180

def build_parser(parser):
    parser.add_argument('clin_flagged', 
                        help='A required input file')
    parser.add_argument('genotype_variants',
                        help='Output file of variants to annotate with varscan readcounts')

def format_variants(series):
    """
    Varscan requires specially formated variant list where
    deletion start position also needs to be [start]-1
    """
    index_cols = ['chrom','start']
    if '-' in series['Var_Base']:
        return pd.Series([series['chrom'],
                          series['start']-1],
                          index=index_cols)
    else:
        return pd.Series([series['chrom'],
                         series['start']],
                         index=index_cols)
 
def action(args):
    
    #Make dataframe of clinically flagged positions 
    cols=['chrom','start','stop','Ref_Base','Var_Base','Clinically_Flagged']
    flagged_variants=pd.read_csv(args.clin_flagged, delimiter='\t', header=None, names = cols, index_col=False)
    #set chr to be a string. 
    flagged_variants['chrom'] = flagged_variants['chrom'].astype('str')

    #input file for varscan, lists postions to call variants at
    genotype_positions = args.genotype_variants

    #Process the clinically flagged positions, format for varscan readcounts function    
    varscan_format_variants = flagged_variants.apply(format_variants, axis = 1)
    clean_df = varscan_format_variants.drop_duplicates(['chrom','start'])

    var_cals = ['chrom','start']
    clean_df.to_csv(genotype_positions, index=False, columns = var_cals, header=None,sep='\t')

