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
from munging.annotation import chromosomes,GenomeIntervalTree, UCSCTable


log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument('refgene', type=Opener(), 
                        help='RefGene file')
    parser.add_argument('bd_file', type=Opener(), 
                        help='Breakdancer output')
    parser.add_argument('-o', '--outfile', type=Opener('w'), metavar='FILE',
                        default=sys.stdout, help='output file')



def set_gene_event(start_pos, chrm, genes):
    """Given start position and gene interval tree, 
    return gene and event"""
    start=int(start_pos)
    matching_genes=genes[chrm].search(start)
    
    if len(matching_genes)<1:
        gene=['Intergenic',]
    else:
        gene=[]
        for start, stop, data in matching_genes:
            gene.append(data['name2'])
    event=chrm+':'+start_pos
    genes=';'.join(str(x) for x in set(gene))

    return event, genes

def action(args):
    genes = GenomeIntervalTree.from_table(args.refgene, parser=UCSCTable.REF_GENE, mode='tx')

    # read in the entire input file so that we can sort it
    in_fieldnames=['#Chr1','Pos1','Orientation1','Chr2','Pos2','Orientation2','Type','Size','Score','num_Reads','num_Reads_lib','Allele_frequency','SampleID']
    reader = csv.DictReader(filter(lambda row: row[0]!='#', args.bd_file), delimiter='\t', fieldnames=in_fieldnames)

    rows = list(reader)
    
    output = []
    for row in rows:
        # each segment is assigned to a gene or exon if either the
        try:
            chr1 = 'chr'+str(chromosomes[row['#Chr1']])
            chr2 = 'chr'+str(chromosomes[row['Chr2']])
        except KeyError:
            print('chrm not being processed: {} or {}'.format(row['#Chr1'], row['Chr2']))
            continue

        row['Event_1'], row['Gene_1']=set_gene_event(row['Pos1'], chr1, genes)
        row['Event_2'], row['Gene_2']=set_gene_event(row['Pos2'], chr2, genes)

        #discard those between -101 and 101
        if int(row['Size']) not in range(-101,101):
            if row['Type']=='CTX':
                row['Size']='N/A'
            output.append(row)

    sorted_output = sorted(output, key=itemgetter('num_Reads'), reverse=True)
    fieldnames=['Event_1','Event_2','Type','Size','Gene_1','Gene_2','num_Reads']
    writer = csv.DictWriter(args.outfile, extrasaction='ignore',fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    writer.writerows(sorted_output)

