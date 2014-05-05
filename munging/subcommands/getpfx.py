"""
Get prefixes files (PFX.[12].fastq.gz/GATKfinal.bam/CNV_bins.txt) for running pipeline.

Usage: 
    munge getpfx /path/to/data

"""

import logging
import os
import glob
import sys

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('datadir', help='Path to directory containing fastq files.')
    parser.add_argument('-s','--separator', help = "separator for list of prefixes",
                        default = '\n')
    
def action(args):
    infiles = glob.glob(os.path.join(args.datadir, '*.gz'))
    prefixes = set(os.path.basename(f).split('.')[0] for f in infiles)
    if len(infiles)<1:
        infiles = glob.glob(os.path.join(args.datadir, '*CNV_bins.txt'))
        prefixes = set(os.path.basename(f).split('_')[0] for f in infiles)
        if len(infiles)<1:
            infiles = glob.glob(os.path.join(args.datadir, '*GATKfinal.bam'))
            prefixes = set(os.path.basename(f).split('.')[0] for f in infiles)
            if len(infiles)<1:
                print 'Need either gz, CNV_bins or GATKfinal.bam to get prefix from'
                sys.exit()

    print args.separator.join(sorted(prefixes))
    
