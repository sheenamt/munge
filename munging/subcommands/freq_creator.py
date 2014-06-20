"""
Calculate tallies of variants and write anovar output

Usage:

    munge freq_creator dbname machine assay -o outfile


"""

import logging
import sqlite3
import csv
import argparse
import sys
from itertools import count, groupby

from munging.utils import dict_factory

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('dbname',
                        help='Name of an sqlite database file')
    parser.add_argument('machine', choices = ['miseq','hiseq'],
                        help = 'name of machine to represent in the output')
    parser.add_argument('assay', choices = ['oncoplex','coloseq'],
                        help = 'name of assay to represent in the output')
    parser.add_argument('-o','--outfile', type = argparse.FileType('w'),
                        default = sys.stdout,
                        help='Name of the output file')

def action(args):

    con = sqlite3.connect(args.dbname)
    # queries return dicts
    con.row_factory = dict_factory

    cur = con.cursor()

    controls = ['LMG-098', 'LMG-098A', 'LMG-098B', 'C066N', 'LMG-240', 'LMG-240A', 'LMG-240B', 'OPX-240', 'OPX-240A', 'OPX-240B','LMG-241', 'CON-0001T', 'CON-0003T', 'LMG-240-1-10', 'CON-0228T', 'CON-762R', 'CON0228T', 'CON762R']
    controlstr = ','.join("'%s'"%c for c in controls)


    # retrieve the variants
    cmd = """
    select *, count(*) as tally from (
    select run, chromosome, start, end, ref_base, var_base
    from variants
    join run_info using (run)
    where pfx not in (%s)
    and machine = ?
    and assay = ?
    group by pfx, chromosome, start, end, ref_base, var_base
    order by run)
    group by chromosome, start, end, ref_base, var_base
    """ % controlstr

    cur.execute(cmd, (args.machine, args.assay,))
    tallies = cur.fetchall()

    cmd = """
    select pfx
    from variants
    join run_info using (run)
    where pfx not in (%s)
    and machine = ?
    and assay = ?
    group by pfx""" % controlstr

    cur.execute(cmd, (args.machine, args.assay,))
    sample_count = len(cur.fetchall())

    writer = csv.DictWriter(
        args.outfile,
        fieldnames = ['chromosome','start','end','ref_base','var_base', 'freq', 'counts'],
        extrasaction = 'ignore',
        delimiter='\t')

    for row in tallies:
        row['counts'] = '[%s/%s]' % (row['tally'], sample_count)
        row['freq'] = '%.4f' % (float(row['tally'])/sample_count)
        writer.writerow(row)
