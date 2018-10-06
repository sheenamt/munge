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
from munging.annotate import (assign, partition, read_refgene,
                               chromosomes, check_overlapping)


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


def get_exons(starts, ends, strand='+'):
    """Return a list of (exon, start, end). `exon` is an int, and `starts`
    and `ends` are strings representing comma-delimited lists of
    chromosomal coordinates. `strand` corresponds to the `strand`
    column of the refGene annotation table; if '-', the order of the
    exon numbers is reversed.

    """

    startpos = [int(i) for i in starts.split(',') if i.strip()]
    endpos = [int(i) for i in ends.split(',') if i.strip()]
    assert len(startpos) == len(endpos)

    numbers = range(1, len(startpos) + 1)

    if strand == '-':
        numbers.reverse()

    # numbers = ['exon{}'.format(n) for n in numbers]

    return zip(numbers, startpos, endpos)


def build_trees(fname):
    """Returns (genes,exone). Each is a dictionary keyed by chromosome
    (`genes`) or gene (`exons`). Values are trees generated by
    `partition`.

    """

    sort_key = itemgetter(1, 2)  # sorts by (start, end)
    genes = defaultdict(list)
    exon_trees = {}
    for row in read_refgene(fname):
        # this fails with a KeyError if a nonstandard chromosome name
        # is encountered (assume these have been filtered out already)

        chrm = str(chromosomes[row['chrom']])
        # gene name from refseq file 
        gene = '{name2}'.format(**row)

        genes[chrm].append((gene, int(row['txStart']), int(row['txEnd'])))

        exons = get_exons(row['exonStarts'], row['exonEnds'], row['strand'])
        exons = sorted(exons, key=sort_key)
        check_overlapping(exons)
        exon_trees[gene] = partition(exons)

        # TODO: if we need to add more features (ie, introns), do it here

    # create a gene tree for each chromosome
    gene_trees = {}
    for chrm, rows in genes.iteritems():
        rows = sorted(rows, key=sort_key)
        check_overlapping(rows)
        gene_trees[chrm] = partition(rows)

    return gene_trees, exon_trees


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
        start1=int(row['Pos1'])
        start2=int(row['Pos2'])
        if str(chr1) in ignored_chrms: # or str(chr2) in ignored_chrms:
            continue
        gene1 = assign(genes[chr1], start1)
        gene2 = assign(genes[chr2], start2)
        row['Event_1'] = 'chr'+row['#Chr1']+':'+row['Pos1']
        row['Event_2'] = 'chr'+row['Chr2']+':'+row['Pos2']

        if gene1:
            row['Gene_1'] = gene1
        else:
            row['Gene_1'] = 'Intergenic'

        if gene2:
            row['Gene_2'] = gene2
        else:
            row['Gene_2'] = 'Intergenic'

        output.append(row)

    fieldnames=['Event_1','Event_2','Type','Size','Gene_1','Gene_2','num_Reads']
    writer = csv.DictWriter(args.outfile, extrasaction='ignore',fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    writer.writerows(output)
