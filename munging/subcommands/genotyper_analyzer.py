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
    parser.add_argument('readcount_output',
                        help='Output of varscan readcounts function')
    parser.add_argument('output', type=argparse.FileType('w'),
                        help='path and name of output to write',
                        default=sys.stdout)
    parser.add_argument('--dec_file',
                        help='annovar formatted file with DEC scores')

def format_indels(series):
    """
    Varscan requires specially formated variant list where
    deletion == DEL-(length of deletion)-Variants deleted
    deletion start position also needs to be [start]-1
    insertion == INS-(length of insertion)-Variants inserted
    """
    if series['start'] == series['stop']:
        position = 'chr'+series['chrom']+':'+str(series['start']) 
    else:
        position = 'chr'+series['chrom']+':'+str(series['start'])+'-'+str(series['stop'])
    index_cols = ['Position','chrom','start','Ref_Base','Var_Base','varscan_start','varscan_ref', 'varscan_variant','Clinically_Flagged']
    if '-' in series['Var_Base']:
        return pd.Series([position,
                          series['chrom'],
                          series['start'],
                          series['Ref_Base'],
                          series['Var_Base'],
                          series['start']-1,
                          '-',
                          '-'.join(['DEL',str(len(series['Ref_Base'])),series['Ref_Base']]),
                          series['Clinically_Flagged']],
                          index=index_cols)
    elif '-' in series['Ref_Base']:
        return pd.Series([position,
                         series['chrom'],
                         series['start'],
                         series['Ref_Base'],
                         series['Var_Base'],
                         series['start'],
                         '-',
                         '-'.join(['INS',str(len(series['Var_Base'])),series['Var_Base']]),
                         series['Clinically_Flagged']],
                         index=index_cols)

    else:
        return pd.Series([position,
                          series['chrom'],
                          series['start'],
                          series['Ref_Base'],
                          series['Var_Base'],
                          series['start'],
                          series['Ref_Base'],
                          series['Var_Base'],
                          series['Clinically_Flagged']],
                         index=index_cols)
 
def parse_varscan_line(line):
    """
    Parse each varscan line into 3 namedtupeles:
    Position
    Reference
    Variant 
    """
    Variant = namedtuple('Variant', ['base','reads','strands','avg_qual','map_qual','plus_reads','minus_reads'])
    Position = namedtuple('Position', ['chrom','position','ref_base', 'depth'])
    Reference = namedtuple('Reference', ['base','reads','strands','avg_qual','map_qual','plus_reads','minus_reads', 'misc'])

    #Split the line by tabs
    sp = line.strip('\n').split('\t')
    # the first three entries are position info, 4th is the depth for minimum quality asked for at readcounts program call
    position = Position(sp[0], sp[1], sp[2], sp[4])
    # column 5 is the reference info, has an extra useless column, is colon delimited
    ref_call = Reference(*sp[5].split(':'))
    # the rest of the columns have more variants, in the same format as the ref info without the extra useless column , each colon delimited
    #If there is only a reference call, sp[6] will be empty string
    if len(sp[6]) > 6:
        #If no specific variant was asked for, all variants have colon delimited info, split with tabs
        if ':' in sp[6]:
            candidates = [s.split(':') for s in sp[6:]]
        #If a specific variant was asked for, it will be columns 6:13, tab delimited
        else:
            query_variant = Variant(*sp[6:13])
            candidates = [s.split(':') for s in sp[13:]]
            candidates.append(query_variant)
            #print candidates
        variants =[Variant(*var) for var in candidates if len(var) > 6]
    else:
        variants = []
    return (position, {'variants': variants, 'reference': ref_call})

    
def action(args):
    
    #Make dataframe of clinically flagged positions 
    cols=['chrom','start','stop','Ref_Base','Var_Base','Clinically_Flagged']
    flagged_variants=pd.read_csv(args.clin_flagged, delimiter='\t', header=None, names = cols, index_col=False)
    #set chr to be a string. 
    flagged_variants['chrom'] = flagged_variants['chrom'].astype('str')

    #output file from varscan readcounts 
    genotype_calls=args.readcount_output

    #final analysis file
    genotype_analysis=args.output

    #Process the clinically flagged positions, format for varscan readcounts function    
    varscan_format_variants = flagged_variants.apply(format_indels, axis = 1)

    # if DEC file provided, add column for DEC values
    if args.dec_file:
        dec_df = pd.read_csv(args.dec_file, sep='\t', header=None, names=['generic', 'uw_dec_pvalue', 'chrom', 'start', 'stop', 'Ref_Base', 'Var_Base'])
        dec_df['chrom'] = dec_df['chrom'].astype('str')
        
        varscan_format_variants = pd.merge(varscan_format_variants, dec_df, how='left', on=['chrom', 'start', 'Ref_Base', 'Var_Base'])
        varscan_format_variants['uw_dec_pvalue'] = varscan_format_variants['uw_dec_pvalue'].fillna(value=-1)

    #parse the varscan output into preferred format 
    reader = open(genotype_calls, 'rU')
    reader.next()

    #For each line in the genotype output, print the info about the position/queried variant
    # and match to the clinically flagged comment 
    for line in reader:
        info = parse_varscan_line(line)
        #info[0] is Position, info[1] is 
        chrom = info[0][0]
        pos_start = int(info[0][1])
        depth = info[0][3]
        varscan_format_variants.loc[(varscan_format_variants['chrom'] == chrom) & (varscan_format_variants['varscan_start'] == pos_start), 'Valid_Reads'] = depth
        varscan_format_variants.loc[(varscan_format_variants['chrom'] == chrom) & (varscan_format_variants['varscan_start'] == pos_start), 'Reference_Reads'] = info[1]['reference'][1]
        for variant in info[1]['variants']:
            varscan_format_variants.loc[(varscan_format_variants['chrom'] == chrom) & (varscan_format_variants['varscan_start'] == pos_start) & (varscan_format_variants['varscan_variant'] == variant[0]), 'Variant_Reads']=variant[1]
        
    #Add column for VAF
    varscan_format_variants['Variant_Read_Fraction'] = varscan_format_variants['Variant_Reads'].astype(float) / varscan_format_variants['Valid_Reads'].astype(float)
    varscan_format_variants['Variant_Read_Fraction'] = varscan_format_variants['Variant_Read_Fraction'].round(4)
    
    header = ['Position','Ref_Base','Var_Base','Clinically_Flagged','Valid_Reads','Reference_Reads', 'Variant_Reads', 'Variant_Read_Fraction']
    if args.dec_file:
        header.append('uw_dec_pvalue')

    varscan_format_variants.to_csv(genotype_analysis, na_rep=0, index=False, columns=header, sep='\t')


    
    
    
