"""
Munges AnnotSV annotation of GRIDSS output
"""
import sys
import subprocess
import logging
import os
import argparse
import pandas as pd

from munging.annotation import multi_split

log = logging.getLogger(__name__)
pd.options.display.width = 180


def build_parser(parser):
    parser.add_argument('annotsv', 
                        help='A required input file')
    parser.add_argument('-o', '--outfile',
                        help='Output file', default=sys.stdout,
                        type=argparse.FileType('w'))


                                                                                
def parse_sv(data):
    ''' Split the ALT and Info fields
    '''
    #Create Event1
    data['Event1']='chr'+str(data['SV chrom'])+':'+str(data['SV start'])

    a,b=multi_split(data['ALT'],'[]')

    #Create Event2
    #the position has a : in it while the sequence does not
    if ':' in a:
        data['Event2']='chr'+a
        data['Seq']=b
    else:
        data['Event2']='chr'+b
        data['Seq']=a

    #Parse read depth and SVtype
    info=dict(item.split('=') for item in data['INFO'].split(";") if "=" in item)
    data['SV_Type']=info['SVTYPE']
    data['Ref_Reads']=info['REFPAIR']
    data['Var_Reads']=info['VF']
    data['EventID']=info['EVENT']
    
    return pd.Series(data)

    
def action(args):
    #Make dataframe of annotsv annotation
    annotsv_df=pd.read_csv(args.annotsv, delimiter='\t', index_col=False, usecols=['SV chrom','SV start','SV end', 'ALT','Gene name','NM','QUAL',
                                                                                   'FILTER','INFO','location','promoters','1000g_event', '1000g_AF', 
                                                                                   'Repeats_coord_left', 'Repeats_type_left', 
                                                                                   'Repeats_coord_right', 'Repeats_type_right',
                                                                                   'Mim Number', 'Phenotypes', 'Inheritance'])

    #Parse the parts we care about
    annotsv_df=annotsv_df.apply(parse_sv, axis=1)

    var_cals = ['Event1', 'Event2', 'SV_Type','Ref_Reads', 'Var_Reads', 'EventID','Gene name','NM','QUAL','FILTER','location','promoters','1000g_event', '1000g_AF', 'Repeats_coord_left', 'Repeats_type_left', 'Repeats_coord_right', 'Repeats_type_right','Mim Number', 'Phenotypes', 'Inheritance']
    annotsv_df.to_csv(args.outfile, index=False, columns=var_cals,sep='\t')

    # combine chr:start-end
    # remove annotations
    # rename annotations
    
