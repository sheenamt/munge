"""Annotate Pindel output with genes, exons, and other features.

The `annotations` file is in the same format as the refGene table, but
has been filtered to contain no overlapping features (ie, using
filter_refseq). Using an unfiltered refGene file will cause an error!

"""

import sys
import argparse
import csv
from collections import defaultdict
from operator import itemgetter
import logging

from munging.utils import Opener
from munging.annotation import (assign,build_trees)


log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument('refgene', type=Opener(), 
                        help='RefGene file, filtered by preferred transcripts file')
    parser.add_argument('pindel_vcfs', action='append', nargs='+',
                        help='Input files which are vcfs from pindel output')
    parser.add_argument('-o', '--outfile', type=Opener('w'), metavar='FILE',
                        default=sys.stdout, help='output file')


def parse_event(data):
    ''' Return the length and type of event '''
    
    #Parse read depth and SVtype
    info=dict(item.split('=') for item in data['INFO'].split(";") if "=" in item)
    size=abs(int(info['SVLEN']))
    if info['SVTYPE'] == 'RPL':
        svtype='DEL'
    else:
        svtype=info['SVTYPE']
    end=info['END']
    return size,info['SVTYPE'], end

def action(args):
    genes,exons = build_trees(args.refgene)
    output = []

    #Skip the header lines 
    (pindel_vcfs,) = args.pindel_vcfs
    for vcf in pindel_vcfs:
        with open(vcf, 'rU')  as f:
            # read in the entire input file so that we can sort it
            fieldnames=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','READS']
            reader = csv.DictReader(filter(lambda row: row[0]!='#', f), delimiter='\t',fieldnames=fieldnames)

            rows = list(reader)
            for row in rows:
                row['Size'], row['Event_Type'],row['End']=parse_event(row)
                if row['Size']<10:
                    continue
                # each segment is assigned to a gene if either the
                # start or end coordinate falls within the feature boundaries.
                chr1=str(row['CHROM'])
                start=int(row['POS'])

                gene1 = assign(genes[chr1], start)
                region = assign(exons[gene1], start)

                out_fieldnames=['Gene','Gene Region','Event_Type','Size','Position','Reads','Transcripts']
        
                row['Position']='chr'+str(chr1)+':'+str(row['POS'])+'-'+str(row['End'])
                row['Reads']=int(row['READS'].split(',')[-1])

                row['Gene Region']=region

                if gene1:
                    row['Gene'] = gene1
                else:
                    row['Gene'] = 'Intergenic'
                output.append(row)

    sorted_output = sorted(output, key=itemgetter('Reads'), reverse=True)  #Sort on reads

    out_fieldnames=['Gene','Gene Region','Event_Type','Size','Position','Reads']
    writer = csv.DictWriter(args.outfile, extrasaction='ignore',fieldnames=out_fieldnames, delimiter='\t')
    writer.writeheader()
    writer.writerows(sorted_output) 
