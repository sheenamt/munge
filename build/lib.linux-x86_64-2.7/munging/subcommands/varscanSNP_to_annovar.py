"""
Clean varscanSNP file, including converting from IUPAC, removing 'chr' prefix, creating read counts

Usage:

munge varscanSNP_to_annovar $PFX.varscanSNP $PFX.varscanSNP.annA
"""
import argparse
import csv
import sys

def build_parser(parser):
    parser.add_argument(
        'infile', type=argparse.FileType('rU'),
        default=sys.stdin)
    parser.add_argument(
        '-o','--outfile',
        help='Output file', default = sys.stdout,
        type = argparse.FileType('w'))


def IUPAC_converter(reader):
    """
    Convert IUPAC to DNA bases, unless converted base is the same as the reference base
    """
    new_var = {'A' : ('A'),
               'C' : ('C'),
               'G' : ('G'),
               'T' : ('T'),
               'R' : ('A','G'),
               'Y' : ('C','T'),
               'M' : ('A','C'),
               'K' : ('G','T'),
               'S' : ('C','G'),
               'W' : ('A','T'),
               'B' : ('C','G','T'),
               'D' : ('A','G','T'),
               'H' : ('A','C','T'),
               'V' : ('A','C','G'),
               }

    for row in reader:
       #skip the header row if present
        if row[0].startswith('Chrom'):
            continue
        chrm=row[0].strip('chr')
        start=row[1]
        stop=row[1]
        ref_base=row[2]
        var_base=row[3]
        reads=row[4]+'|'+row[5]
        var_qual=row[10]
        #for each var_base, convert using the new_var table
        # unless base is N, in which case, skip it. 
        if var_base=='N':
            continue
        for b in new_var[var_base]:
            #don't list row if the new var_base is the same as ref_base
            if not b is ref_base:
                new_row=chrm, start, stop, ref_base, b, reads, var_qual
            else:
                continue
            yield new_row

def action(args):
    fname = args.infile
    outname= args.outfile
    reader = csv.reader(fname, delimiter="\t")
    writer = csv.writer(outname,
                        delimiter='\t',
                        )
    
    for row in IUPAC_converter(reader):
        writer.writerow(row)
 
