"""
Crawl analysis files to create one analysis file with all info

Usage:

 munge combined_msi $SAVEPATH -o $SAVEPATH/${RUN}_Combined_MSI

"""

import logging
import os
import csv
import sys
import argparse
import collections
from itertools import count, groupby, chain, ifilter

from munging.utils import dict_factory
from munging.annotation import split_chr_loc
from munging import filters
from munging.utils import walker, munge_pfx

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('path',
                        help='Path to analysis files')
    parser.add_argument('-o','--outfile', type = argparse.FileType('w'),
                        default = sys.stdout,
                        help='Name of the output file')


def tally_msi(data, file_name):
    """Count the total MSI sites found and
    count the mutated sites
    """
    total, mutants = 0, 0
    with open(os.path.join(file_name), 'rU') as fname:
        reader = csv.DictReader(fname, delimiter='\t')
        pfx = file_name.split('.')[0].split('/')[-1].strip('_msi')
        reader = sorted(reader, key=itemgetter('Position'))
        for row in data:
            for key, group in groupby(reader, key=itemgetter('Position')):
                #for each position in control (row[0])
                if row['Position'] == key:
                    for info in group:
                        if int(info['Avg_read_depth']) >= 30:
                            value = float(row['Ave']) + (2 * float(row['Std']))
                            if int(info['Number_Peaks']) >= value:
                                total += 1
                                mutants += 1
                            else:
                                total += 1

    return total, mutants, pfx
    
def action(args):

    files = ifilter(filters.any_analysis, walker(args.path))
    files = ifilter(filters.msi_analysis, files)
    #sort the files so that the output in the workbook is sorted
    files=sorted(files)

    variant_keys = ['Position', 'Gene']
    specimens, annotation, prefixes=parse_pindel(variant_keys, files, args.path)
    annotation_headers = [
        'Gene_Region',
        'Event_Type',
        'Size',
        'Transcripts',
        ]

    writer = csv.DictWriter(args.outfile, fieldnames = variant_keys + annotation_headers + prefixes,  extrasaction = 'ignore', delimiter = '\t')
    writer.writeheader()
    for variant in sorted(specimens.keys()):
        d = {k:v for k,v in zip(variant_keys,variant)}
        d.update({pfx:specimens[variant].get(pfx) for pfx in prefixes})
        d.update(annotation[variant])
        writer.writerow(d)
