"""
Parse variant files from pipeline and NIST to create QC Variant file

Usage:

munge qc_variant pipeline_file NIST_file -o output file 
"""
import argparse
import csv
import sys

def build_parser(parser):
    parser.add_argument(
        'pipeline', type=argparse.FileType('rU'), 
        default=sys.stdin)
    parser.add_argument(
        'NISTfile', type=argparse.FileType('rU'),
        default=sys.stdin)
    parser.add_argument(
        '-o','--outfile',
        help='Output file', default = sys.stdout,
        type = argparse.FileType('w'))


#read in bed file: chr, start, stop
#parse annA file: 	#if this position is captured by a particular probe
#if chr is the same:
    #if ann[start] >= bed[start] AND ann[stop] <= bed[stop]:
       #print entire ann row to outfile 


def match(pipe, nist):
    """
    Return list of variants found in all three input files (NIST, LMG/OPX-240 output
    Matches if chrm, start, stop, ref_base, and var_base are the same.
    """
    matchedlist = []
    for nist_line in nist:
        for pipe_line in pipe:    
            try:
                if nist_line[3:8]==pipe_line[3:8]:
                    if nist_line[1:8] not in matchedlist:
                        matchedlist.append(nist_line[1:8])
            except IndexError:
                print 'Index Error', nist_line, pipe_line
                continue
    return matchedlist

def action(args):

    pipe = list(csv.reader(args.pipeline, delimiter="\t"))
    nist = list(csv.reader(args.NISTfile, delimiter="\t"))

    output = match(pipe, nist)
    
    writer = csv.writer(args.outfile,
                        csv.QUOTE_NONE,
                        delimiter='\t')

    for i in output:
        writer.writerow(i)
    





