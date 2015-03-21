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


from munging.utils import dict_factory
from munging.annotation import fix_pfx

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

    # retrieve the variants
    cmd = """
    select  chromosome, start, end, ref_base, var_base, count(distinct pfx) as tally
    from variants
    join run_info using (run, project)
    where not is_control
    and machine = ?
    and assay = ?
    group by chromosome, start, end, ref_base, var_base
    """

    cur.execute(cmd, (args.machine, args.assay,))
    tallies = cur.fetchall()


    cmd = """
    select count(distinct pfx) as sample_count
    from variants
    join run_info using (run, project)
    where not is_control
    and machine = ?
    and assay = ?
    """

    cur.execute(cmd, (args.machine, args.assay,))
    sample_count = cur.fetchone()['sample_count']

    writer = csv.DictWriter(
        args.outfile,
        fieldnames = ['chromosome','start','end','ref_base','var_base', 'freq', 'counts'],
        extrasaction = 'ignore',
        delimiter='\t')

    for row in tallies:
        row['counts'] = '[%s/%s]' % (row['tally'], sample_count)
        row['freq'] = '%.4f' % (float(row['tally'])/sample_count)
        writer.writerow(row)
