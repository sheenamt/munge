
"""
Filter a file of genomic positions given ranges of start positions
specified in a bed file

Usage:

munge bed_parser bedfile input_file -o output_file

"""
from collections import defaultdict
from operator import itemgetter
import csv
import sys

def build_parser(parser):
    parser.add_argument(
        'bedfile',
        help='Path to Bed File, used as filter criteria',
        default=sys.stdin)
    parser.add_argument(
        'annfile',
        help='Path to file that needs to be filtered ',
        default=sys.stdin)
    parser.add_argument(
        '-o','--outfile',
        help='Output file',
        default = sys.stdout)


def coords(row, chromosome=0, start=1, stop=2):
    """Return chromosome, start, end - assuming that column indices are
    the same in bedfile and giantfile.

    """

    return row[chromosome], int(row[start]), int(row[stop])

def action(args):

    bedfile, giantfile, outfile = args.bedfile, args.annfile, args.outfile
    # prepare a dictionary of chromosome: set(positions)
    # includes all positions between start-stop so we
    # can test for membership efficiently
    ranges = defaultdict(set)
    with open(bedfile) as f:
        for row in csv.reader(f, delimiter='\t'):
            row[0]=row[0].strip('chr')
            chr, beg, end = coords(row)
            ranges[chr].update(range(beg, end + 1))
    # now we can filter
    with open(giantfile) as g, open(outfile,'w') as o:
        writer = csv.writer(o, delimiter='\t')
        # is the start position among the allowed posisions for this
        # chromosome?
        for row in csv.reader(g, delimiter='\t'):
            chr, beg, _ = coords(row)
            if beg in ranges[chr]:
                writer.writerow(row)

