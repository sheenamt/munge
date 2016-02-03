"""
Parse Illumina Amplicon output to create amplicon quality file

Usage:

munge amplicon_metrics amplicon_bed fastq_dir outdir --top_output $PROJECT.Combined_Amplicons.txt
"""
import argparse
import sys
import csv
import pandas as pd
import glob
from itertools import ifilter
from os import path

from munging import parsers,filters
from munging.utils import walker

def build_parser(parser):
    parser.add_argument(
        'amplicons', type=argparse.FileType('rU'),
        default=sys.stdin)
    parser.add_argument(
        'run_metrics_dir', help='run directory')
    parser.add_argument(
        'project',
        default=sys.stdin)
    parser.add_argument(
        '--top_output', type=argparse.FileType('w'),
        default=sys.stdout)

def action(args):
    #get run-amplicon file
    files = ifilter(filters.amplicon_coverage, walker(args.run_metrics_dir))
    files = sorted(files)
    assert len(files) is 1
    for pth in files:
        with open(path.join(pth.dir, pth.fname)) as fname:
            run_metrics=pd.read_csv(fname, delimiter='\t')
            #the target column header is empty, so add a name
            run_metrics.rename(columns={'Unnamed: 0':'Target'}, inplace=True)
            #clean up the metrics file to remove empty columns
            clean_metrics = run_metrics.dropna(axis='columns',how='all')

            #Grab samples from file, removing 'Target' columns 
            samples=list(clean_metrics.columns.values)
            samples.remove('Target')

    #grab amplicon bed
    amplicons=pd.read_csv(args.amplicons)

    #merge metrics on the amplicon targets
    merged = pd.merge(amplicons,clean_metrics, how='inner', left_on='Target', right_on='Target')

    #Print top-level-output
    merged.to_csv(args.top_output, index=False,sep='\t')

    #Print sample level output
    for sample in samples:
        pfx=sample+'-'+args.project
        header=['Target','Gene','Position']
        header.append(sample)
        sample_info=merged[header]
        #Expected : project/pfx/pfx.Amplicon_Analysis.txt
        sample_out = path.join('output',pfx,pfx+'.Amplicon_Analysis.txt')
        sample_info.to_csv(sample_out,index=False, sep='\t')
