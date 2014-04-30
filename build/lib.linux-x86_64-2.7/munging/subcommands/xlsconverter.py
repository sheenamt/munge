"""
Create xls workbook from all output files

usage:

 munge xlsmaker (Analysis/Combined) $SAVEPATH/$PFX* -o $SAVEPATH/${PFX}_Analysis.xls

"""
import csv
import argparse
import sys
import re
import os

from xlwt import Workbook


book = Workbook()

def build_parser(parser):
    parser.add_argument('type', choices = ['Analysis','Combined'],
        help = 'Type of files to create workbook, Analysis or Combined')
    parser.add_argument(
        'infiles', action = 'append', nargs = '+',
        help='Input files')
    parser.add_argument(
        '-o','--outfile',
        help='Output file', default = sys.stdout,
        type = argparse.FileType('w'))


def float_if_possible(strg):
    """
    Convert integers to float instead of string where applicable. 
    """
    try:
        return float(strg)
    except ValueError:
        return strg

def action(args):
    filetype = args.type 
    print filetype
    (infiles, ) = args.infiles
    for fname in infiles:
        (f_path, f_name) = os.path.split(fname)
        if re.search(str(filetype), f_name):
            (f_short_name, f_extension) = os.path.splitext(f_name)
            if len(f_short_name) >= 30:
                f_short_name=f_short_name[:30]
            print f_short_name
            sheet= book.add_sheet(f_short_name)
            Reader = csv.reader(open(fname, 'rU'), delimiter='\t')
            for rowx, row in enumerate(Reader):
                for colx, value in enumerate(row):
                    sheet.write(rowx, colx, float_if_possible(value))

    book.save(args.outfile)



