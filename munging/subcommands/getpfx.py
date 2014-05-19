"""
Get prefixes files (PFX.[12].fastq.gz) for running pipeline.

Usage: 
    munge getpfx /path/to/data

"""

import logging
import os
import glob

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('datadir', help='Path to directory containing fastq files.')
    parser.add_argument('-s','--separator', help = "separator for list of prefixes",
                        default = '\n')
    
def action(args):
    infiles = glob.glob(os.path.join(args.datadir, '*.gz'))
    prefixes = set(os.path.basename(f).split('.')[0] for f in infiles)
    print args.separator.join(sorted(prefixes))
    
