"""
Parse variant files from pipeline, 1000G, and Complete Genomics to create QC Variant file

Usage:

munge qc_variant pipeline_file kg_file Complete_genomics_file -o output file 
"""
import argparse
import csv
import sys

def build_parser(parser):
    parser.add_argument(
        'pipeline', type=argparse.FileType('rU'), 
        default=sys.stdin)
    parser.add_argument(
        'kgfile', type=argparse.FileType('rU'),
        default=sys.stdin)
    parser.add_argument(
        'CGfile', type=argparse.FileType('rU'),
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


def match(pipe, kg, cg):
    """
    Return list of variants found in all three input files (1000G, Complete Genomes, LMG/OPX-240 output
    Matches if chrm, start, stop, ref_base, and var_base are the same.
    """
    matchedlist = []
    for cg_line in cg:
        for kg_line in kg:    
            for pipe_line in pipe:    
                try:
                    if (cg_line[3:8]==kg_line[3:8]) and (kg_line[3:8]==pipe_line[3:8]):
                        if cg_line[1:8] not in matchedlist:
                            matchedlist.append(cg_line[1:8])
                
                except IndexError:
                    print 'Index Error', cg_line, kg_line, pipe_line
                    continue
    return matchedlist

def action(args):

    pipe = list(csv.reader(args.pipeline, delimiter="\t"))
    kg = list(csv.reader(args.kgfile, delimiter="\t"))
    cg = list(csv.reader(args.CGfile, delimiter="\t"))

    output = match(pipe, kg, cg)
    
    writer = csv.writer(args.outfile,
                        csv.QUOTE_NONE,
                        delimiter='\t')

    for i in output:
        writer.writerow(i)
    





