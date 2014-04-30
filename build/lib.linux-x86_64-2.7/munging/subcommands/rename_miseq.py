"""
Rename MiSeq files for pipeline processing

Usage:
     munge rename_miseq /path/to/data /dest/path/ 

"""

import logging
import shutil
import os
import glob
import re

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('datadir', 
                        help='Source directory')
    parser.add_argument('-d','--outdir', 
                        help='Destination directory [default %(default)s]',
                        default = '.')
    
def action(args):

    infiles = glob.glob(os.path.join(args.datadir, '*.gz'))
    for f in infiles:
        prefix = os.path.basename(f).split('_')[0]
        if re.search('R1', f):
            shutil.move(f, '%s.1.fastq.gz' % (prefix))
        elif re.search('R2', f):
            shutil.move(f, '%s.2.fastq.gz' % (prefix))

