"""Annotate Pindel output with genes, exons, and other features.

The `annotations` file is in the same format as the refGene table, but
has been filtered to contain no overlapping features (ie, using
filter_refseq). Using an unfiltered refGene file will cause an error!

"""

import sys
import argparse
import logging
import pandas as pd
from munging.utils import Opener
from munging.annotation import chromosomes,GenomeIntervalTree, UCSCTable, multi_split, merge_transcript_data

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument('refgene', 
                        help='RefGene file, filtered by preferred transcripts file')
    parser.add_argument('pindel_vcfs', action='append', nargs='+',
                        help='Input files which are vcfs from pindel output')
    parser.add_argument('--multi_reads', action='store_true',
                        help='Expect bbmerged and bwamem reads')
    parser.add_argument('-o', '--outfile', type=Opener('w'), metavar='FILE',
                        default=sys.stdout, help='output file')


def parse_vcf_info(info):
    ''' Return the length and type of event, corrects end position if necessary '''
    #Parse read depth and SVtype
    info_dict = dict(item.split('=') for item in info.split(";") if "=" in item)
    # report absolute value for size
    size = abs(int(info_dict['SVLEN']))
    if info_dict['SVTYPE'] == 'RPL':
        svtype = 'DEL'
    else:
        svtype = info_dict['SVTYPE']
    # To maintain consistency with the old Perl annotations, add 1 to the end
    # This means insertions will be listed with consecutive breakends (e.g. 16:68771418-68771419)
    # Deletions and tandem duplications will be listed with breakends separated by a number
    # of positions equal to the size of the event (e.g. 16:817016-817032 for a size 15 deletion)
    # In all cases, the breakends represent the UNAFFECTED nucleotides closest to the variation
    end = int(info_dict['END']) + 1
    return pd.Series([size, svtype, end])

def get_annotations(row, genome_tree):
    chrom = chromosomes[row['CHROM']]
    chrom_tree = genome_tree[chrom]
    start = int(row['POS'])
    end = int(row['End'])
    start_transcripts = [ x[2] for x in chrom_tree[start] ]
    end_transcripts = [ x[2] for x in chrom_tree[end] ]
    spanning_transcripts = [ x[2] for x in chrom_tree[start:end + 1] ]

    start_genes, _, start_annotations = merge_transcript_data(start_transcripts, start, report_utr=True)
    end_genes, _, end_annotations = merge_transcript_data(end_transcripts, end, report_utr=True)
    _, region_types , _ = merge_transcript_data(spanning_transcripts, start, end + 1, report_utr=True)

    genes = sorted(set(start_genes + end_genes))
    if len(genes) == 0:
        genes = ['Intergenic']
    annotations = sorted(set(start_annotations + end_annotations))

    if 'EXON' in region_types:
        gene_region = 'EXONIC'
    elif 'UTR' in region_types:
        gene_region = 'UTR'
    elif 'INTRON' in region_types:
        gene_region = 'INTRONIC'
    else:
        gene_region = 'Intergenic'

    gene_label = ';'.join(str(x) for x in genes)
    transcript_label = ';'.join(str(x) for x in annotations) #transcripts #combine_transcripts(set(transcripts)) #

    return pd.Series([gene_label, transcript_label, gene_region])

def parse_vcf_reads(read_info):
    return int(read_info.split(',')[-1])

def get_position(row):
    chrom = chromosomes[row['CHROM']]
    start = str(row['POS'])
    end = str(row['End'])
    return "{}:{}-{}".format(chrom, start, end)

# def ranges(i):
#     for a, b in itertools.groupby(enumerate(i), lambda (x, y): y - x):
#         b = list(b)
#         yield b[0][1], b[-1][1]
    
def action(args):
    #Create interval tree of Transcripts, grouped by chr
    gt = GenomeIntervalTree.from_table(args.refgene)

    # specify columns to import and names
    if args.multi_reads:
        headers=['CHROM','POS','INFO','bbmergedREADS','bwamemREADS']
        cols_to_use = [0,1,7,9,10]
    else:
        headers=['CHROM','POS','INFO','READS']
        cols_to_use = [0,1,7,9]

    # import VCFs into DataFrames
    (pindel_vcfs,) = args.pindel_vcfs
    readers = []
    for vcf in pindel_vcfs:
        with open(vcf, 'rU')  as f:
            reader = pd.read_csv(f, comment='#', delimiter='\t', header=None, usecols=cols_to_use, names=headers)
        readers.append(reader)

    # concatenate DataFrames into one
    df = pd.concat(readers, ignore_index=True)
    chroms = chromosomes.keys()
    # filter out entries containing events on unsupported chromosomes
    df = df[df['CHROM'].isin(chroms)]
    # parse the contents of the INFO field
    df[['Size', 'Event_Type', 'End']] = df['INFO'].apply(parse_vcf_info)
    # do not include small insertion/deletion calls from Pindel
    df = df[df['Size'] > 10]
    # add annotations
    df[['Gene', 'Transcripts', 'Gene_Region']] = df.apply(get_annotations, genome_tree=gt, axis=1)
    # create the Position field
    df['Position'] = df.apply(get_position, axis=1)

    # parse the contents of the reads field(s)
    if args.multi_reads:
        df['bbmergedReads'] = df['bbmergedREADS'].apply(parse_vcf_reads)
        df['bwamemReads'] = df['bwamemREADS'].apply(parse_vcf_reads)
    else:
        df['Reads']=df['READS'].apply(parse_vcf_reads)

    # sort and save to file
    if args.multi_reads:
        df = df.sort_values(['bwamemReads', 'Position'], ascending=[False, True])  #Sort on reads
        out_fields=['Gene','Gene_Region','Event_Type','Size','Position','bbmergedReads', 'bwamemReads','Transcripts']
    else:
        df = df.sort_values(['Reads', 'Position'], ascending=[False, True])  #Sort on reads
        out_fields=['Gene','Gene_Region','Event_Type','Size','Position','Reads', 'Transcripts']

    df.to_csv(args.outfile, sep='\t', index=False, columns=out_fields)

