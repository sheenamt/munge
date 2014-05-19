"""
Crawl control files to find std and ave of each position across control run

Usage:

munge msi_pipeline_samples  /path/to/control/file /path/to/sample/files -o output_file

"""

import os
import csv
import sys
import argparse
from itertools import groupby
from operator import itemgetter


def build_parser(parser):
    parser.add_argument('control_file', type=argparse.FileType('rU'),
                        default=sys.stdin,
                        help='Path to control file')
    parser.add_argument('sample_files', nargs='+',
                        help='Path to sample files')
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'),
                        default=sys.stdout,
                        help='Name of the output file')


def tally_msi(data, file_name):
    """Count the total MSI sites found and
    count the mutated sites
    """
    total, mutants = 0, 0
    with open(os.path.join(file_name), 'rU') as fname:
        reader = csv.DictReader(fname, delimiter='\t')
        pfx = file_name.split('.')[0].split('/')[-1]
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


def write_output(writer, total, mutants, pfx):
    """Write formatted output to file
    """
    totals = ["Total microsatellite loci: %s " % total]
    muts = ["Total unstable loci: %s " % mutants]
    if total:
        freq = (float(mutants) / total)
    else:
        freq=0
    mut_freq = ["Fraction of unstable loci for %s: %0.4f " % (pfx, freq)]

    writer.writerow(totals)
    writer.writerow(muts)
    writer.writerow(mut_freq)


def action(args):
    control_info = csv.DictReader(args.control_file, delimiter='\t')
    writer = csv.writer(args.outfile, delimiter='\t')
    #Store the dictreader in a variable to loop through it twice
    data = [row for row in control_info]
    for file_name in args.sample_files:
        total, mutants, pfx = tally_msi(data, file_name)
        write_output(writer, total, mutants, pfx)
