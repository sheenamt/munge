"""
Parse bedtools output for intervals not meeting given threshold, annotate with RefSeq data

Usage:

munge coverage_metrics /path/to/bedtools/output /path/to/refseq -o /path/to/output
"""
import sys
import csv
import pandas as pd
import munging.annotation as ann
from munging.utils import Opener
from operator import itemgetter
from natsort import natsorted

def build_parser(parser):
    parser.add_argument('target_coverage',
                        help='Bedtools output for base level coverage')
    parser.add_argument('refgene', 
                        help='RefGene file')
    parser.add_argument('-t','--threshold', type=int, default=50,
                        help='Threshold for outputing failed coverage')
    parser.add_argument('-o', '--outfile', type=Opener('w'), metavar='FILE',
                        default=sys.stdout, help='output file')


def action(args):
    #read in refgene into genome interval tree
    gt = ann.GenomeIntervalTree.from_table(args.refgene)

    output = []

    #read in bedfile
    headers=['CHROM','START','1_BASED_POSITION','READS']
    data_types={'CHROM': str, 'START': int, '1_BASED_POSITION' : int, 'READS' : int}
    df = pd.read_csv(args.target_coverage, comment='#', delimiter='\t',header=None, usecols=[0,1,4,5],names=headers, dtype=data_types)

    df['POS']=df['START'] + (df['1_BASED_POSITION']-1)

    #remove all rows that pass the thresholds
    failing_bases = df[df['READS']<=args.threshold]
    
    #Match the chrm annotation of the exons file
    #Now, groupy by START, mean READS, and label interval with first and last POS
    data=failing_bases.groupby(['CHROM','START'],as_index=False).agg({'READS':['mean','median'],
                                                                        'POS':['min','max']})

    #Convert to a dictionary for processing clearly
    rows = data.to_dict(orient='records')
    
    # add annotations
    for row in rows:
        # determine genomic position
        chrom=ann.chromosomes[row[('CHROM','')]]
        start=int(row[('POS','min')])
        end=int(row[('POS','max')])
        
        # set the position label
        if start==end:
            row['Position'] = '{}:{}'.format(chrom, start)
        else:
            row['Position'] = '{}:{}-{}'.format(chrom, start, end)
            
        # find all transcripts in the interval [start, end]
        transcripts = [ x[2] for x in gt[chrom][start:end + 1] ]

        # create annotations
        genes = ann.gene_info_from_transcripts(transcripts)
        region_types = ann.region_info_from_transcripts(transcripts, start, end + 1, report_utr=True)
        transcript_labels = ann.transcript_info_from_transcripts(transcripts, start, end + 1, report_utr=True)

        row['Mean Coverage'] = round(row[('READS','mean')], 1)
        row['Median Coverage'] = round(row[('READS','median')], 1)
        row['Gene'] =';'.join(genes)
        row['Gene_Region']=';'.join(region_types)
        row['Transcripts']=';'.join(transcript_labels)
        output.append(row)
        
    output=natsorted(output,key=itemgetter('Position'))  # sort on position 
    out_fieldnames=['Position','Gene','Gene_Region','Mean Coverage', 'Median Coverage','Transcripts']
    writer = csv.DictWriter(args.outfile, extrasaction='ignore',fieldnames=out_fieldnames, delimiter='\t')
    writer.writeheader()
    writer.writerows(output) 
                
