"""
Parse Illumina Amplicon output to create amplicon quality file

Usage:

munge probe_coverage perprobecoverage -o output
"""
import argparse
import sys
import csv
import pandas as pd
from munging.utils import Opener
import pandas
from operator import itemgetter

def build_parser(parser):
    parser.add_argument('perprobecoverage', type=argparse.FileType('rU'),
                        help='path to probe coverage from bedtools ')
    parser.add_argument('-o', '--outfile', type=Opener('w'), metavar='FILE',
                        default=sys.stdout, help='output file')

def action(args):

    output = []
    headers=['chrom','start','stop', 'Probe', 'strand', 'mean']
    reader = pandas.read_csv(args.perprobecoverage, comment='#', delimiter='\t',header=None,usecols=[0,1,2,3,5], names=headers)
    
    rows = reader.T.to_dict().values()
        
    for row in rows:
        row['Position']='{}:{}-{}'.format(row['chrom'],row['start'],row['stop'])
        row['MeanCoverage']=int(round(row['mean']))
        output.append(row)

    output.sort(key=itemgetter('Probe'))
    out_fieldnames=['Position','Probe','MeanCoverage']

    writer = csv.DictWriter(args.outfile, extrasaction='ignore',fieldnames=out_fieldnames, delimiter='\t')
    writer.writeheader()
    writer.writerows(output) 

                                              
