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
    ignored_chrms = args.ignore_chrms or []

    # read in the entire input file so that we can sort it
    in_fieldnames=['#Chr1','Pos1','Orientation1','Chr2','Pos2','Orientation2','Type','Size','Score','num_Reads','num_Reads_lib','Allele_frequency','SampleID']
    reader = csv.DictReader(filter(lambda row: row[0]!='#', args.bd_file), delimiter='\t', fieldnames=in_fieldnames)

    rows = list(reader)

    output = []
    for row in rows:
        # each segment is assigned to a gene or exon if either the
        # start or end coordinate falls within the feature boundaries.
        chr1=str(row['#Chr1'])
        chr2=str(row['Chr2'])

        if str(chr1) in ignored_chrms or str(chr2) in ignored_chrms:
            continue

        start1=int(row['Pos1'])

        try:
            gene1 = assign(genes[chr1], start1)
            region1= assign(exons[gene1], start1)
        except KeyError:
            gene1='Intergenic or off target'
            region1='Intergenic or off target'
        row['Event_1'] = 'chr'+chr1+':'+row['Pos1']
        row['Gene_1_Region']=region1

        if gene1:
            row['Gene_1'] = gene1
        else:
            row['Gene_1'] = 'Intergenic'


        start2=int(row['Pos2'])

        try:
            gene2 = assign(genes[chr2], start2)
            region2= assign(exons[gene1], start1)
        except KeyError:
            gene='Intergenic or off target'
            region2='Intergenic or off target'

        row['Event_2'] = 'chr'+chr2+':'+row['Pos2']
        row['Gene_2_Region']=region2
        if gene2:
            row['Gene_2'] = gene2
        else:
            row['Gene_2'] = 'Intergenic'
        #discard those between -101 and 101
        if int(row['Size']) not in range(-101,101):
            output.append(row)

    fieldnames=['Event_1','Event_2','Type','Size','Gene_1','Gene_1_Region','Gene_2','Gene_2_Region','num_Reads']
    writer = csv.DictWriter(args.outfile, extrasaction='ignore',fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    writer.writerows(output)
