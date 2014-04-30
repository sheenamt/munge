"""
Crawl control files to find std and ave of each position across control run

Usage:

munge msi_control /path/to/control/files -o output_file

"""

import os
import csv 
import sys
import argparse 
from collections import defaultdict
from itertools import groupby, ifilter
from operator import itemgetter
from numpy import std, array, average

from munging import filters
from munging.utils import walker



def build_parser(parser):
    parser.add_argument('path', 
                        help='Path to analysis files')
    parser.add_argument('-o','--outfile', type = argparse.FileType('w'),
                        default = sys.stdout,
                        help='Name of the output file')

def action(args):
    control=defaultdict(list)
    # apply a series of filters to files
    files = ifilter(filters.any_analysis, walker(args.path))
    #sort the files so that the output in the workbook is sorted
    files=sorted(files)
    count=0
    for pth in files:
        with open(os.path.join(args.path, pth.fname)) as fname:
            reader = csv.DictReader(fname, delimiter='\t')
            reader=sorted(reader, key=itemgetter('Position'))
            for k, g in groupby(reader, key=itemgetter('Position')):
                for row in g:
                    if int(row['Avg_read_depth']) >= 30:
                        control[row['Position']].append(int(row['Number_Peaks']))
    header=['Position', 'Std', 'Ave', 'Count']
    writer = csv.writer(args.outfile, quoting=csv.QUOTE_MINIMAL, delimiter='\t')
    writer.writerow(header)
    for k,v in sorted(control.items()):
        count=len(v)
        a=array(v)
        std=a.std()
        ave=average(a)
        row=[k, std, ave, count]
        writer.writerow(row)
