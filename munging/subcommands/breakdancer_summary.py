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
import pandas

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
    event=str(chrm)+':'+str(start_pos)
    genes=';'.join(str(x) for x in set(gene))

    return event, genes

def action(args):
    genes = GenomeIntervalTree.from_table(args.refgene, parser=UCSCTable.REF_GENE, mode='tx')

    #try to match chr string in reference file and input data
    chrm=False
    if 'chr' in genes.keys()[0]:
        chrm = True

    # read in only the columns we care about, because real data can be too large sometimes
    headers=['#Chr1','Pos1','Chr2','Pos2','Type','Size','num_Reads']
    reader = pandas.read_csv(args.bd_file, comment='#', delimiter='\t',header=None,usecols=[0,1,3,4,6,7,9], names=headers)
    #Convert to a dictionary for processing clearly
    rows = reader.T.to_dict().values()
    output = []
    for row in rows:
        # each segment is assigned to a gene or exon if either the
        #only normal chr are process, GL amd MT are ignored
        try:
            if chrm:
                chr1 = 'chr'+str(chromosomes[row['#Chr1']])
                chr2 = 'chr'+str(chromosomes[row['Chr2']])
            else:
                chr1 = str(chromosomes[row['#Chr1']])
                chr2 = str(chromosomes[row['Chr2']])
        except KeyError:
            print('chrm not being processed: {} or {}'.format(row['#Chr1'], row['Chr2']))
            continue

        row['Event_1'], row['Gene_1']=set_gene_event(row['Pos1'], chr1, genes)
        row['Event_2'], row['Gene_2']=set_gene_event(row['Pos2'], chr2, genes)

        #discard those between -101 and 101
        if int(row['Size']) not in range(-101,101):
            #display absolute value for size
            row['Size'] = abs(int(row['Size']))
            
            if row['Type']=='CTX':
                row['Size']='N/A'

            output.append(row)

    output.sort(key=itemgetter('Event_1'))
    output.sort(key=itemgetter('num_Reads'), reverse=True)
    fieldnames=['Event_1','Event_2','Type','Size','Gene_1','Gene_2','num_Reads']
    writer = csv.DictWriter(args.outfile, extrasaction='ignore',fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    writer.writerows(output)

