"""
Description: convert UCSC refgene.txt files to BED format 
"""
 
import os
import sys 
 
def build_parser(parser):
    parser.add_argument('--input', nargs='?', default=sys.stdin)
    parser.add_argument('--output', nargs='?', help='Output file', default = sys.stdout)

 
def get_int_list(l):
    return [int(i) for i in l.strip(',').split(',')]
 
def get_string_list(a):
    return ','.join([str(i) for i in a])
 
def bed_key(d):
    return([d['chrom'], int(d['chromStart']), int(d['chromEnd'])])

def action(args):
    genes = {}
    for line in open(args.input, 'r'):
        ls = line.strip().split('\t')
        starts, stops =  get_int_list(ls[9]), get_int_list(ls[10])
        lengths = get_string_list([stop-start for stop,start in zip(stops, starts)])
        relstarts = get_string_list([start - int(ls[4]) for start in starts])
        # For format see http://genome.ucsc.edu/FAQ/FAQformat.html#format1
        features = ['chrom','chromStart','chromEnd','name', 'refseq','score','strand','thickStart',
                    'thickEnd','itemRgb','exonCount','exonSizes','exonStarts'] 
        gene_entry = dict([('chrom', ls[2].strip('chr')),
                           ('chromStart', ls[4]),
                           ('chromEnd', ls[5]),
                           ('name', ls[12]),
                           ('refseq', ls[1]),
                           ('score', 1),
                           ('strand', ls[3]),
                           ('thickStart', ls[4]),
                           ('thickEnd', ls[5]),
                           ('itemRgb', 0),
                           ('exonCount', ls[8]),
                           ('exonSizes', lengths),
                           ('exonStarts', relstarts)])
        refseq = ls[1]
        
        # Ensure that each refseq is only in the table once
        if (refseq not in genes 
            and 'NM' in refseq):
            genes[refseq] = gene_entry

    output=open(args.output,'w')
    for gene in sorted(genes.values(), key=bed_key):
        output.write('\t'.join([str(gene[f]) for f in features]) + '\n')
    output.close()
