"""Annotate BreakDancer output with genes, exons, and other features.

The `annotations` file is in the same format as the refGene table, but
has been filtered to contain no overlapping features (ie, using
filter_refseq). Using an unfiltered refGene file will cause an error!

"""

import sys
import argparse
import logging
from munging.utils import Opener
import munging.annotation as ann
import pandas as pd

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument('refgene', type=Opener(), 
                        help='RefGene file')
    parser.add_argument('bd_file', type=Opener(), 
                        help='Breakdancer output')
    parser.add_argument('-g', '--genes', type=Opener(),
                        help='Text file of genes to flag')
    parser.add_argument('-o', '--outfile', type=Opener('w'), metavar='FILE',
                        default=sys.stdout, help='output file')

def add_genes(chr_pos_pair, genome_tree):
    """returns a string with the concatenated set of genes found at the chr_pos_pair coordinates in genome_tree"""
    chrom = ann.chromosomes[chr_pos_pair[0]]
    pos = int(chr_pos_pair[1])
    transcripts = [x[2] for x in genome_tree[chrom][pos]]
    genes = ann.gene_info_from_transcripts(transcripts)
    return ';'.join(genes)

def add_event(row):
    chrom = str(row[0])
    pos = str(row[1])
    return "{}:{}".format(chrom,pos)

def filter_genes(row, gene_set):
    """returns True if any gene in the row is found in the gene_set; otherwise returns False"""
    genes1 = row['Gene_1'].split(';')
    genes2 = row['Gene_2'].split(';')
    row_genes = set(genes1 + genes2)
    gene_intersection = row_genes & gene_set
    return len(gene_intersection) > 0

def action(args):
    gt = ann.GenomeIntervalTree.from_table(args.refgene)

    # read in only the columns we care about, because real data can be too large sometimes
    headers=['Chr1','Pos1','Chr2','Pos2','Type','Size','num_Reads']
    df = pd.read_csv(args.bd_file, comment='#', delimiter='\t',header=None, usecols=[0,1,3,4,6,7,9], names=headers)

    # discard rows with size in [-101, 101]
    df = df[df['Size'].abs() > 101]
    # update Size when Type is CTX
    df.loc[df['Type'] == 'CTX', 'Size'] = 'N/A'

    # discard rows containing events on unsupported chromosomes
    chroms = ann.chromosomes.keys()
    df = df[df['Chr1'].isin(chroms) & df['Chr2'].isin(chroms)]

    # check that any calls remain before applying functions or pandas with crash
    if len(df) > 0:

        # add events columns
        df['Event_1'] = df[['Chr1', 'Pos1']].apply(add_event, axis=1)
        df['Event_2'] = df[['Chr2', 'Pos2']].apply(add_event, axis=1)
        # add genes columns
        df['Gene_1'] = df[['Chr1', 'Pos1']].apply(add_genes, genome_tree = gt, axis=1)
        df['Gene_2'] = df[['Chr2', 'Pos2']].apply(add_genes, genome_tree = gt, axis=1)

        if args.genes:
            # read in genes to keep
            gene_df = pd.read_csv(args.genes, comment='#', delimiter='\t',header=None, usecols=[0], names=['gene'])
            gene_set = set(gene_df['gene'])
            # discard rows that don't contain a gene from the gene_set
            df['Flagged_Genes'] = df.apply(filter_genes, gene_set=gene_set, axis=1)

        # sort and save the relevant columns
        df.sort_values(['num_Reads','Event_1'], ascending=[False, True], inplace=True)

    if args.genes:
        out_fields = ['Event_1','Event_2','Type','Size','Gene_1','Gene_2', 'Flagged_Genes', 'num_Reads']
    
    else:
        out_fields = ['Event_1','Event_2','Type','Size','Gene_1','Gene_2', 'num_Reads']
    
    df.to_csv(args.outfile, index=False, sep='\t', columns=out_fields)

