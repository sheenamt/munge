
"""
Parse annovar input files to contain only regions captured in BED file

Usage:

munge annovar_bed_parser bedfile $SAVEPATH/$PFX.annA -o $SAVEPATH/$PFX.ann

"""
import argparse
import csv
import sys

def build_parser(parser):
    parser.add_argument(
        'bedfile', type=argparse.FileType('rU'), 
        default=sys.stdin)
    parser.add_argument(
        'annfile', type=argparse.FileType('rU'),
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


def match(bedinfo, anninfo):
    """
    Return list of variants found in both the bed file and the annA file in pipeline. 
    Matches if chr is the same, ann[start] >= bed[start] AND ann[stop] <= bed[stop]
    """
    matchedlist = []
    for annline in anninfo:
        for bedline in bedinfo:    
            bedline[0]=bedline[0].strip('chr')   
            try:
                if (bedline[0]==annline[0]):
                    if (int(annline[1])>=int(bedline[1])) and (int(annline[1])<=int(bedline[2])):
                        if annline not in matchedlist:
                            matchedlist.append(annline)
            except IndexError:
                print 'Index Error', bedline, annline
                continue
    return matchedlist

def action(args):

    bedinfo = list(csv.reader(args.bedfile, delimiter="\t"))
    anninfo = list(csv.reader(args.annfile, delimiter="\t"))
    
    output = match(bedinfo, anninfo)
    
    writer = csv.writer(args.outfile,
                        csv.QUOTE_NONE,
                    delimiter='\t')
    for i in output:
        writer.writerow(i)
    





