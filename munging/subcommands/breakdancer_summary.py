"""Annotate BreakDancer output with genes, exons, and other features.

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
from munging.annotation import (assign, partition, read_refgene,
                                chromosomes, check_overlapping,
                                get_exons, build_trees)


log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument('refgene', type=Opener(), 
                        help='RefGene file, filtered by preferred transcripts file')
    parser.add_argument('bd_file', type=Opener(), 
                        help='Breakdancer output')
    parser.add_argument('-i','--ignore-chrms', nargs='+',
                        help='CHRMs to ignore for CNV annotation')
    parser.add_argument('-o', '--outfile', type=Opener('w'), metavar='FILE',
                        default=sys.stdout, help='output file')




def action(args):
    genes, exons = build_trees(args.refgene)

    # read in the entire input file so that we can sort it
    in_fieldnames=['#Chr1','Pos1','Orientation1','Chr2','Pos2','Orientation2','Type','Size','Score','num_Reads','num_Reads_lib','Allele_frequency','SampleID']
    reader = csv.DictReader(filter(lambda row: row[0]!='#', args.bd_file), delimiter='\t', fieldnames=in_fieldnames)

    rows = list(reader)

    output = []
    for row in rows:
        # each segment is assigned to a gene or exon if either the
        try:
            chr1 = str(chromosomes[row['#Chr1']])
            chr2 = str(chromosomes[row['Chr2']])
        except KeyError:
            print('chrm not being processed: {} or {}'.format(row['#Chr1'], row['Chr2']))
            continue

        start1=int(row['Pos1'])
        try:
            gene1 = assign(genes[chr1], start1)
            region1 = assign(exons[gene1], start1)
        except KeyError:
            gene1='Intergenic'
            region1='Intergenic'
        row['Event_1'] = 'chr'+chr1+':'+row['Pos1']
        row['Gene_1'] = gene1
        row['Gene_1_Region']=region1

        start2=int(row['Pos2'])
        try:
            gene2 = assign(genes[chr2], start2)
            region2= assign(exons[gene2], start2)
        except KeyError:
            gene2='Intergenic'
            region2='Intergenic'
        row['Event_2'] = 'chr'+chr2+':'+row['Pos2']
        row['Gene_2'] = gene2
        row['Gene_2_Region']=region2

        #discard those between -101 and 101
        if int(row['Size']) not in range(-101,101):
            if row['Type']=='CTX':
                row['Size']='N/A'
            output.append(row)

    sorted_output = sorted(output, key=itemgetter('num_Reads'), reverse=True)
    fieldnames=['Event_1','Event_2','Type','Size','Gene_1','Gene_1_Region','Gene_2','Gene_2_Region','num_Reads']
    writer = csv.DictWriter(args.outfile, extrasaction='ignore',fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    writer.writerows(sorted_output)
