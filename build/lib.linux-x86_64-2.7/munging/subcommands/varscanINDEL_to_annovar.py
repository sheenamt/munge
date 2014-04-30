"""
Clean varscan file, including converting indel notation to bases, removing 'chr' prefix, creating read counts

Usage:

munge varscanINDEL_to_annovar $PFX.varscan $PFX.varscanINDEL.annA
"""
import argparse
import csv
import sys
from munging.annotation import multi_split

def build_parser(parser):
    parser.add_argument(
        'infile', type=argparse.FileType('rU'),
        default=sys.stdin)
    parser.add_argument(
        '-o','--outfile',
        help='Output file', default = sys.stdout,
        type = argparse.FileType('w'))


def indel_cleaner(reader):
    """
    Remove the unwanted info in the variant base column of the varscan file for processing through Annovar
    """
    for row in reader:
       #skip the header row if present
        if row[0].startswith('Chrom'):
            continue
        chrm=row[0].strip('chr')
        start=row[1]
        stop=row[1]
        ref_base=row[2]
        #split the var_base to remove the unwanted info
        var_base=multi_split(row[3], '*/')
        reads=row[4]+'|'+row[5]
        var_qual=row[10]
        #create the new row with just the bases in the var_base 
        if var_base[0].startswith('-'):
            ref_base=var_base[0].strip('-')
            #calculate the stop position of the deletion
            #this is how Steve has is set up
            stop=int(start)+len(ref_base)
            start=int(start)+1            
            var_base='-'
        elif var_base[0].startswith('+'):
            var_base=var_base[0].strip('+')
            ref_base='-'
        new_row=chrm, start, stop, ref_base, var_base, reads, var_qual
        yield new_row

def action(args):
    fname = args.infile
    outname= args.outfile
    reader = csv.reader(fname, delimiter="\t")
    writer = csv.writer(outname,
                        delimiter='\t',
                        )
    for row in indel_cleaner(reader):
        writer.writerow(row)
 
