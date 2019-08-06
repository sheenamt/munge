"""
Parse bedtools output for intervals not meeting given threshold, annotate with RefSeq data

Usage:

munge coverage_metrics /path/to/bedtools/output /path/to/refseq -o /path/to/output
"""
import sys
import csv
import pandas as pd
from munging.utils import Opener
from munging.annotation import chromosomes,GenomeIntervalTree, UCSCTable, define_transcripts

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
    #read in refgene into interval tree
    #Create interval tree of introns and exons,  grouped by chr
    exons = GenomeIntervalTree.from_table(open(args.refgene, 'r'), parser=UCSCTable.REF_GENE, mode='exons')


    #try to match chr string in reference file and input data
    chrm=False
    if 'chr' in exons.keys()[0]:
        chrm = True

    output = []

    #read in bedfile
    with open(args.target_coverage, 'rU')  as f:
        headers=['CHROM','START','END','GENE','1_BASED_POSITION','READS']
        df = pd.read_csv(f, comment='#', delimiter='\t',header=None,usecols=[0,1,4,5],names=headers)
        df['POS']=df['START']+(df['1_BASED_POSITION']-1)

        #remove all rows that pass the thresholds
        failing_bases = df[df['READS']<=args.threshold]
        
        
        #Match the chrm annotation of the exons file
        #Now, groupy by START, mean READS, and label interval with first and last POS
        data=failing_bases.groupby(['CHROM','START'],as_index=False).agg({'READS':['mean','median'],
                                                                          'POS':['min','max']})
        if chrm:
            data.loc[:,'CHROM']='chr'+data['CHROM'].astype(str)
        else:
            data.loc[:,'CHROM']=data['CHROM'].astype(str)

        #NOW, ANNOTATE:

        #Convert to a dictionary for processing clearly
        rows = data.T.to_dict().values()
        
        for row in rows:
            chrm=row[('CHROM','')]
            start=int(row[('POS','min')])
            end=int(row[('POS','max')])
            
            #This is a snp
            if start==end:
                chrm_exons=exons[chrm].search(start)
                row['Position'] = '{}:{}'.format(chrm, start)
            else:
                chrm_exons=exons[chrm].search(start,end)
                row['Position'] = '{}:{}-{}'.format(chrm, start, end)
                
            #now we either have an interval or we don't
            if chrm_exons:
                gene1, region, transcripts = define_transcripts(chrm_exons)


            else:
                gene1 = ['Intergenic']
                region = ['Intergenic']
                transcripts = []

            row['Mean Coverage'] = row[('READS','mean')]
            row['Median Coverage'] = row[('READS','median')]
            row['Gene'] =';'.join(str(x) for x in set(gene1))
            row['Gene_Region']=';'.join(str(x) for x in set(region))
            row['Transcripts']=';'.join(str(x) for x in set(transcripts))
            output.append(row)
            
        out_fieldnames=['Position','Gene','Gene_Region','Mean Coverage', 'Median Coverage','Transcripts']
        writer = csv.DictWriter(args.outfile, extrasaction='ignore',fieldnames=out_fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(output) 
                
