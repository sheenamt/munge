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
    parser.add_argument('specimen_folder',
                        help='A required specimen ID folder')
    parser.add_argument('varscan',
                        help='Absolute path to varscan jar file')
    parser.add_argument('output', type=argparse.FileType('w'),
                        help='path and name of output to write',
                        default=sys.stdout)
    parser.add_argument('--min_coverage', default=8, type=int,
                        help='Minimum read depth at a position to make a call')
    parser.add_argument('--min_base_qual', default=30, type=int,
                        help='Minimum read depth at a position to make a call')
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


 
def run_varscan_readcounts(mpileup, varscan, genotype_positions, genotype_calls, min_qual, min_coverage):
    """ 
    Run the varscan readcounts function
    """

    cmd = ['java', '-jar', varscan, 'readcounts', mpileup, '--variants-file', genotype_positions, '--min-base-qual', str(min_qual), '--min-coverage', str(min_coverage), '--output-file', genotype_calls]
    varscan_readcounts = subprocess.Popen(cmd)
    varscan_readcounts.wait()

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
    # the first three entries are position info
    position = Position(sp[0], sp[1], sp[2], sp[4])
    # column 5 is the reference info, has an extra useless column
    ref_call = Reference(*sp[5].split(':'))
    # the rest of the columns have more variants, in the same format as the ref info without the extra useless column 
    candidates = [s.split(':') for s in sp[6:]]
    variants =[Variant(*var) for var in candidates if len(var) > 6]

    return (position, {'variants': variants, 'reference': ref_call})

def action(args):
    
    #Make dataframe of clinically flagged positions 
    cols=['chrom','start','stop','Ref_Base','Var_Base','Clinically_Flagged']
    flagged_variants=pd.read_csv(args.clin_flagged, delimiter='\t', header=None, names = cols, index_col=False)
    #set chr to be a string. 
    flagged_variants['chrom'] = flagged_variants['chrom'].astype('str')

    #setup the sample files. 
    sample_path = args.specimen_folder
    sample_id = os.path.basename(sample_path)

    #mpileup that readcounts runs on 
    mpileup = os.path.join(sample_path, sample_id+'.mpileup')

    #input file for varscan, lists postions to call variants at
    genotype_positions = os.path.join(sample_path,sample_id+'.genotype_input')

    #output file from varscan readcounts 
    genotype_calls=os.path.join(sample_path,sample_id+'.genotype_output') 

    #final analysis file
    genotype_analysis=args.output

    #Process the clinically flagged positions, format for varscan readcounts function    
    varscan_format_variants = flagged_variants.apply(format_indels, axis = 1)

    var_cals = ['chrom','varscan_start','varscan_ref','varscan_variant']
    varscan_format_variants.to_csv(genotype_positions, index=False, columns = var_cals, header=None,sep='\t')

    #Run varscan readcount, creates output file we process next
    run_varscan_readcounts(mpileup, args.varscan, genotype_positions, genotype_calls, args.min_base_qual, args.min_coverage)

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
        varscan_format_variants.loc[(varscan_format_variants['chrom'] == chrom) & (varscan_format_variants['varscan_start'] == pos_start), 'Reference_Reads']=info[1]['reference'][1]
        varscan_format_variants.loc[(varscan_format_variants['chrom'] == chrom) & (varscan_format_variants['varscan_start'] == pos_start), 'Variant_Reads']=0
        for variant in info[1]['variants']:
            print 'variant:', variant
            varscan_format_variants.loc[(varscan_format_variants['chrom'] == chrom) & (varscan_format_variants['varscan_start'] == pos_start) & (varscan_format_variants['varscan_variant'] == variant[0]), 'Variant_Reads']=variant[1]

    header = ['Position','Ref_Base','Var_Base','Clinically_Flagged','Valid_Reads','Reference_Reads','Variant_Reads']
    
    varscan_format_variants.to_csv(genotype_analysis, na_rep='0',index=False,columns=header,sep='\t')


    
    
    
