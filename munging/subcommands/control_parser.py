"""
Compare quality control variants to OPX-240 output to check quality of run 

Usage:

munge control_parser /path/to/HapMap_variant_function /path/to/control/sample_Analysis.txt -o HapMap_QC_Analysis.txt
"""
import argparse
import csv
import sys

import pandas as pd

def build_parser(parser):
    parser.add_argument(
        'control', type=argparse.FileType('rU'), 
        default=sys.stdin)
    parser.add_argument(
        'run_output', type=argparse.FileType('rU'),
        default=sys.stdin)
    parser.add_argument(
        '-o','--outfile',
        help='Output file', default = sys.stdout,
        type = argparse.FileType('w'))


def match(control, run_info):
    """
    Make a list and keep count of variants found in both the qc file and the run output for LMG/OPX-240 sample
    Matches if chr, start are the same. 
    """
    matchedlist = []
    missedlist = []
    concount=0
    for conline in control:
        concount=concount+1
        for runline in run_info:
            if runline[0].startswith('Position'):
                continue
            else:
                pos=multi_split(runline[0], 'chr:-,')
                #if chrm, position, refbase and varbase match
                if (conline[2]==pos[0]) and (int(conline[3])==int(pos[1])) and (conline[5]==runline[1]) and (conline[6] == runline[2]):
                    if conline not in matchedlist:
                        matchedlist.append(conline)
                else:
                   missedlist.append(conline) 
    return matchedlist, missedlist, concount

def action(args):
    #Read in the data from the HapMap VCF file
    cols=['variant_type','anno','chr','start','Ref_Base','Var_Base']
    control_df = pd.read_csv(args.control, delimiter='\t', header=None, usecols=[0,1,2,3,5,6],names = cols, index_col=False)
    control_df['Position']='chr'+control_df['chr'].astype('str')+':'+control_df['start'].astype('str')
    #Get the count of the number of control snps
    expected_count=len(control_df)
    
    #Read in the data from the pipeline data of this sample
    pipeline_df = pd.read_csv(args.run_output, delimiter='\t', header=0, usecols=[0,1,2,4],index_col=False)

    #Merge by finding all keys in the 'left/control' data, whether present in the pipeline data or not
    matched=pd.merge(control_df, pipeline_df, how='left', on=['Position','Ref_Base','Var_Base'], indicator=True)
    
    #Tally the total variants that were found in both outputs
    found=len(matched[matched['_merge']=='both'])

    #Missed data is defined by 'left_only' (ie wasn't found the the 'right/pipeline' data
    missed=matched.loc[(matched['_merge'] == 'left_only')]

    #drop the columns we no longer care about, which makes it easier to add the flavor text for output 
    missed=missed.drop(columns=['chr','start','Variant_Type'])

    output=args.outfile
    headers = ['Position', 'Ref_Base', 'Var_Base', 'variant_type', 'anno']
    #If the 'missed' dataframe is empty, we found all of the variants
    if missed.empty:
        new_line="Found {} of expected {} variants".format(len(matched.index), len(control_df.index))
        output.write(new_line)
    else:
        new_line={'Position':"Found {} of expected {} variants. Missed:".format(found, expected_count),  'Ref_Base':'', 'Var_Base':'', 'Variant_Type':'', 'anno':'', '_merge':''}  #[]#+ [' ']*(len(missed.columns)-1)
        missed.loc[-1]=new_line
        missed.index = missed.index +1
        missed=missed.sort_index()
        missed.to_csv(output, sep='\t', columns=headers, header=None, index=False, mode='a')
