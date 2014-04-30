"""
Compare quality control variants to OPX-240 output to check quality of run 

Usage:

munge control_parser /home/genetics/genetics_pipeline/doc/control_variants /path/to/control/sample_Analysis.txt -o OPX-240_QC_Analysis.txt
"""
import argparse
import csv
import sys

from munging.annotation import multi_split

def build_parser(parser):
    parser.add_argument(
        'control', type=argparse.FileType('rU'), 
        default=sys.stdin)
    parser.add_argument(
        'run_output', type=argparse.FileType('rU'),
        default=sys.stdin)
    parser.add_argument(
        '-o','--outfile',
        help='Output file', default = sys.stdout,
        type = argparse.FileType('w'))


def match(control, run_info):
    """
    Make a list and keep count of variants found in both the qc file and the run output for LMG/OPX-240 sample
    Matches if chr, start are the same. 
    """
    matchedlist = []
    nonmatch= []
    concount=0
    for conline in control:
        concount=concount+1
        for runline in run_info:
            if runline[0].startswith('Position'):
                continue
            else:
                runline=multi_split(runline[0], 'chr:-,')
                if (conline[2]==runline[0]) and (int(conline[3])==int(runline[1])):
                    if conline not in matchedlist:
                        matchedlist.append(conline)
    return matchedlist, concount

def action(args):

    control = list(csv.reader(args.control, delimiter="\t"))
    run_info = list(csv.reader(args.run_output, delimiter="\t"))

    output, concount= match(control, run_info)

    writer = csv.writer(args.outfile,
                        csv.QUOTE_NONE,
                        delimiter='\t', 
                        lineterminator='\n')
    expected = ("Total number of variants expected:", concount)
    actual = ("Total number of variants found:", len(output))
    missed=  ("Variants missed (if any):", )
    writer.writerow(expected)
    writer.writerow(actual)
    writer.writerow(missed)
    for conline in control:
        if conline not in output:
            writer.writerow(conline)
    found=  ("Variants found:", )
    writer.writerow(found)
    for line in output:
        new_line=line[2], line[3], line[1]
        writer.writerow(new_line)
    





