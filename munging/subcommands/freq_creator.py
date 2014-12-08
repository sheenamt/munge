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
    
    controls = ('LMG098','LMG098A','LMG098B','C066N','LMG240','LMG240110'
                'LMG240A','LMG240B','OPX240','OPX240A','OPX240B','LMG241')
    controlstr = ','.join("'%s'"%c for c in controls)


    # retrieve the variants
    cmd = """
    select *, count(*) as tally from (
    select run, chromosome, start, end, ref_base, var_base
    from variants
    join run_info using (run)
    where pfx not in {}
    and pfx not like 'CON%'
    and pfx not like '%NA12878%'
    and machine = ?
    and assay = ?
    group by pfx, chromosome, start, end, ref_base, var_base
    order by run)
    group by chromosome, start, end, ref_base, var_base
    """.format(controls)

    cur.execute(cmd, (args.machine, args.assay,))
    tallies = cur.fetchall()


    cmd = """
    select pfx
    from variants
    join run_info using (run)
    where pfx not in {}
    and pfx not like 'CON%'
    and pfx not like '%NA12878%'
    and machine = ?
    and assay = ?
    group by pfx""".format(controls)

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
