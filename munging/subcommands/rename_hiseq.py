"""
Rename and compress HiSeq files.

Usage:
    munge rename_hiseq /path/to/data /dest/path/ 

"""

import logging
import shutil
import os
import glob

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('datadir', 
                        help='Source directory')
    parser.add_argument('-d','--outdir', 
                        help='Destination directory [default %(default)s]',
                        default = '.')
    parser.add_argument('-n','--dry-run', 
                        help='Dry run (print commands and exit)', dest ='run',
                        action = 'store_false', default = True)
def action(args):
    infiles = glob.glob(os.path.join(args.datadir, '*.txt*'))
    for f in infiles:
        pfx, i, _ = os.path.basename(f).split('_', 2)
        newname = os.path.join(args.outdir, '%(pfx)s.%(i)s.fastq.gz' % locals())
        print newname
        if args.run:
            shutil.move(f, newname)
