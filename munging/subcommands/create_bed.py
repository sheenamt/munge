"""
Parse bed file from agilent into 120bp probes or CONTRA dies

Usage:

munge create_bed /path/to/bed/file -o new_bed_file

"""

import os
import csv 
import sys
import argparse 
from collections import defaultdict
from itertools import groupby, ifilter
from operator import itemgetter

from munging import filters
from munging.utils import walker



def build_parser(parser):
    parser.add_argument('bed', type=argparse.FileType('rU'), 
                        help='Bed file',
                        default=sys.stdin)
    parser.add_argument('-o','--outfile', type = argparse.FileType('w'),
                        default = sys.stdout,
                        help='Name of the output file')

#take in bed file
#format:
#chr# \t start \t stop \t gene 
# for each row:
# return chr# \t start \t start+120 \t gene
#until start >= stop

def action(args):
    bed=args.bed
    reader = csv.reader(bed, delimiter='\t')
    new_probes=[]
    for row in reader:
#        new_probes.append(row)
        start=int(row[1])
        stop=int(row[2])
        while start<stop:
            new_row=(row[0],start,int(start)+120,row[3],'+','+')
            new_probes.append(new_row)
            start=int(start)+40
    writer=csv.writer(args.outfile, delimiter='\t')
    for i in new_probes:
        writer.writerow(i)
