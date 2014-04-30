"""
Parse picard and CNV output to create quality metrics file

Usage:

munge quality_metrics $SAVEPATH/$PFX_CNV_bins.txt $SAVEPATH/$PFX.quality_metrics -o $SAVEPATH/$PFX_quality_metrics_Analysis.txt 
"""
import argparse
import csv
import sys
import re
from numpy import std, array, average

def build_parser(parser):
    parser.add_argument(
        'cnv_file', type=argparse.FileType('rU'), 
        default=sys.stdin)
    parser.add_argument(
        'qm_file', type=argparse.FileType('rU'), 
        default=sys.stdin)
    parser.add_argument(
        '-o','--outfile',
        help='Output file', default = sys.stdout,
        type = argparse.FileType('w'))

def get_values(info):
    """
    Get the standard deviation from the tumor.rd.ori column in the CNV file
    """
    tumors=[]
    for i in info:
        tumors.append(float(i[16])) #this is the tumor.rd.ori column
    b=array(tumors)    
    return b 

def action(args):

    cnv_info = csv.reader(args.cnv_file, delimiter="\t")
    metric_info = csv.reader(args.qm_file, delimiter="\t")
    cnv_info.next()
    writer = csv.writer(args.outfile, quoting=csv.QUOTE_MINIMAL, delimiter='\t')
    a=get_values(cnv_info)
    stdev=["Standard deviation of read depth: %0.2f " %(a).std()]
    ave=["Average read depth: %0.2f " %average(a)]
    writer.writerow(stdev)
    writer.writerow(ave)
    for line in metric_info:
        try:
            if re.search('LIBRARY', line[0]):
                writer.writerow(line)
                writer.writerow(next(metric_info))
        except IndexError:
            continue
            
