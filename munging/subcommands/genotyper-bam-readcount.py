"""
Parses varscan readcounts for specific positions 
"""
import sys
import subprocess
import tempfile
import logging
import shutil
import os
import argparse
import pandas as pd
import csv 
from collections import namedtuple

log = logging.getLogger(__name__)
pd.options.display.width = 180

def build_parser(parser):
    parser.add_argument('bam',
                        help='Sample bam file')
    parser.add_argument('clin_flagged', 
                        help='A required input file')
    parser.add_argument('ref_fasta', 
                        help='reference sequence in the fasta format')
    parser.add_argument('output', type=argparse.FileType('w'),
                        help='path and name of output to write',
                        default=sys.stdout)
    parser.add_argument('--min_mapping_qual', default=0, type=int,
                        help='minimum mapping quality of reads used use the read for counting.')
    parser.add_argument('--min_base_qual', default=0, type=int,
                        help='minimum base quality at a position to use the read for counting.')

def format_indels(series):
    """
    BAM requires specially formated variant list where
    deletion == -Variants
    insertion == +Variants
    """
    if series['start'] == series['stop']:
        position = 'chr'+series['chrom']+':'+str(series['start']) 
    else:
        position = 'chr'+series['chrom']+':'+str(series['start'])+'-'+str(series['stop'])
    index_cols = ['Position','chrom','start','Ref_Base','Var_Base','bam_start','bam_ref', 'readcount_variant','Clinically_Flagged']
    #Deletion
    if '-' in series['Var_Base']:
        return pd.Series([position,
                          series['chrom'],
                          series['start'],
                          series['Ref_Base'],
                          series['Var_Base'],
                          series['start'],
                          series['Ref_Base'],
                          '-'+str(series['Ref_Base']),
                          series['Clinically_Flagged']],
                          index=index_cols)
    #Insertion
    elif '-' in series['Ref_Base']:
        return pd.Series([position,
                          series['chrom'],
                          series['start'],
                          series['Ref_Base'],
                          series['Var_Base'],
                          series['start'],
                          series['Ref_Base'],
                          '+'+str(series['Var_Base']),
                          series['Clinically_Flagged']],
                         index=index_cols)
        
    else:
        #SNP
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

def run_bam_readcounts(bam, clin_flagged, ref_fasta, readcounts_output , min_mapping_qual, min_base_qual):
    """ 
    Run the bam readcounts function
    """
    cmd = ['bam-readcount', '-f', ref_fasta, bam, '-l', clin_flagged, '-q', str(min_mapping_qual), '-b', str(min_base_qual), '--max-warnings', '0']
    f = open(readcounts_output,"wb")
    subprocess.check_call(cmd,stdout=f)

def parse_readcount_line(line):
    """
    Parse each readcount line into 3 namedtupeles:
    Position
    Variant 
    Reference
    """
    Position = namedtuple('Position', ['chrom','position','ref_base', 'depth'])
    Variant = namedtuple('Variant',['base','count','avg_mapping_quality','avg_basequality','avg_se_mapping_quality','num_plus_strand','num_minus_strand','avg_pos_as_fraction','avg_num_mismatches_as_fraction','avg_sum_mismatch_qualities','num_q2_containing_reads','avg_distance_to_q2_start_in_q2_reads','avg_clipped_length','avg_distance_to_effective_3p_end'])
    Reference = namedtuple('Reference', ['base','count','avg_mapping_quality','avg_basequality','avg_se_mapping_quality','num_plus_strand','num_minus_strand','avg_pos_as_fraction','avg_num_mismatches_as_fraction','avg_sum_mismatch_qualities','num_q2_containing_reads','avg_distance_to_q2_start_in_q2_reads','avg_clipped_length','avg_distance_to_effective_3p_end'])

    #Split the line by tabs
    sp = line.strip('\n').split('\t')
    # the first three entries are position info
    position = Position(*sp[:4])
    # the rest of the columns have the ref and  variants
    candidates = [s.split(':') for s in sp[4:]]
    variants =[Variant(*var) for var in candidates if len(var) > 6]
    ref_call = [Reference(*var) for var in candidates if var[0]==sp[2]]
    return (position, {'variants': variants}, {'reference': ref_call})


def action(args):
    
    #Make dataframe of clinically flagged positions 
    cols=['chrom','start','stop','Ref_Base','Var_Base','Clinically_Flagged']
    flagged_variants=pd.read_csv(args.clin_flagged, delimiter='\t', header=None, names = cols, index_col=False)
    #set chr to be a string. 
    flagged_variants['chrom'] = flagged_variants['chrom'].astype('str')

    #bam that readcounts runs on 
    bam = args.bam
    bam_readcounts = args.bam+('.readcounts.txt')

    # #Run varscan readcount, creates output file we process next
    run_bam_readcounts(bam, args.clin_flagged, args.ref_fasta, bam_readcounts,  args.min_mapping_qual, args.min_base_qual)

    #Process the clinically flagged positions, format for varscan readcounts function    
    bamcount_format_variants = flagged_variants.apply(format_indels, axis = 1)
    bam_variants = args.bam+('variants.txt')
    bamcount_format_variants.to_csv(bam_variants, index=False,sep='\t')

    #parse the bam readcount output into preferred format 
    reader = open(bam_readcounts, 'rU')

    #For each line in the genotype output, print the info about the position/queried variant
    # and match to the clinically flagged comment 
    for line in reader:
        info = parse_readcount_line(line)
        #info[0] is Position, info[1] is 
        chrom = info[0][0]
        pos_start = int(info[0][1])

        depth = info[0][3]
        bamcount_format_variants.loc[(bamcount_format_variants['chrom'] == chrom) & (bamcount_format_variants['bam_start'] == pos_start), 'Valid_Reads']= depth
        bamcount_format_variants.loc[(bamcount_format_variants['chrom'] == chrom) & (bamcount_format_variants['bam_start'] == pos_start), 'Reference_Reads']=info[2]['reference'][0][1]
        bamcount_format_variants.loc[(bamcount_format_variants['chrom'] == chrom) & (bamcount_format_variants['bam_start'] == pos_start), 'Variant_Reads']=0
        for variant in info[1]['variants']:
            bamcount_format_variants.loc[(bamcount_format_variants['chrom'] == chrom) & (bamcount_format_variants['bam_start'] == pos_start) & (bamcount_format_variants['readcount_variant'] == variant[0]), 'Variant_Reads']=variant[1]

    header = ['Position','Ref_Base','Var_Base','Clinically_Flagged','Valid_Reads','Reference_Reads','Variant_Reads']
    bamcount_format_variants.to_csv(args.output, index=False,columns=header,sep='\t')


    
    
    
