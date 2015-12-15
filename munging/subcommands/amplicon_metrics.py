"""
Parse Illumina Amplicon output to create amplicon quality file

Usage:

munge data/HHv1.csv amplicon_metrics_file $SAVEPATH/ -o $SAVEPATH/$PFX_amplicon_metrics_Analysis.txt
"""
import argparse
import sys
import csv
import pandas as pd

from os import path

from munging import parsers

def build_parser(parser):
    parser.add_argument(
        'amplicons', type=argparse.FileType('rU'),
        default=sys.stdin)
    parser.add_argument(
        'run_metrics', type=argparse.FileType('rU'),
        default=sys.stdin)
    parser.add_argument(
        '-o', '--outfile',
        help='Output file', default=sys.stdout,
        type=argparse.FileType('w'))

def action(args):
    amplicons=pd.read_csv(args.amplicons)
    run_metrics=pd.read_csv(args.run_metrics, delimiter='\t')

    #the target column header is empty, so add a name
    run_metrics.rename(columns={'Unnamed: 0':'Target'}, inplace=True)

    #clean up the metrics file to remove empty columns
    clean_metrics = run_metrics.dropna(axis='columns',how='all')

    #Grab samples from file, removing 'Target' columns 
    samples=list(clean_metrics.columns.values)
    samples.remove('Target')

    #merge metrics on the amplicon targets
    merged = pd.merge(amplicons,clean_metrics, how='inner', left_on='Target', right_on='Target')
    #Print top-level-output
    merged.to_csv('output.csv', index=False)
    #Print sample level output
    for sample in samples:
        print "sample:", sample
        header=['Target','Gene','Position']
        header.append(sample)
        sample_info=merged[header]
        sample_info.to_csv(sample+'.amplicon_analysis.txt', index=False)
