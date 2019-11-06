"""Annotate Pindel output with genes, exons, and other features.

The `annotations` file is in the same format as the refGene table, but
has been filtered to contain no overlapping features (ie, using
filter_refseq). Using an unfiltered refGene file will cause an error!

"""

import sys
import argparse
import csv
from collections import defaultdict
from operator import itemgetter
import logging
import pandas
from munging.utils import Opener
from munging.annotation import chromosomes,GenomeIntervalTree, UCSCTable, define_transcripts, multi_split
import itertools

csv.field_size_limit(10000000)
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


def parse_event(data):
    ''' Return the length and type of event, corrects end position if necessary '''
    #Parse read depth and SVtype
    info=dict(item.split('=') for item in data['INFO'].split(";") if "=" in item)
    size=int(info['SVLEN'])
    if info['SVTYPE'] == 'RPL':
        svtype='DEL'
    else:
        svtype=info['SVTYPE']
    #Pindel reports insertions wrong by not setting the end position correctly. It may do it with other data, so test based on reported size rather than svtype
    if abs(size)>1 and int(data['POS']) == int(info['END']):
        end=int(info['END'])+1
    else:
        end=int(info['END'])
    return size,svtype, end

def ranges(i):
    for a, b in itertools.groupby(enumerate(i), lambda (x, y): y - x):
        b = list(b)
        yield b[0][1], b[-1][1]
    
def action(args):
    #Create interval tree of introns and exons,  grouped by chr
    exons = GenomeIntervalTree.from_table(open(args.refgene, 'r'), parser=UCSCTable.REF_GENE, mode='exons')
    
    #try to match chr string in reference file and input data
    chrm=False
    if 'chr' in exons.keys()[0]:
        chrm = True

    output = []

    (pindel_vcfs,) = args.pindel_vcfs
    for vcf in pindel_vcfs:
        with open(vcf, 'rU')  as f:
            if args.multi_reads:
                #Skip the header lines, read in only the columns we need because some unnecessary columns can be millions of characters long
                headers=['CHROM','POS','INFO','bbmergedREADS','bwamemREADS']
                reader = pandas.read_csv(f, comment='#', delimiter='\t',header=None,usecols=[0,1,7,9,10], names=headers)
            else:
                #Skip the header lines, read in only the columns we need because some unnecessary columns can be millions of characters long
                headers=['CHROM','POS','INFO','READS']
                reader = pandas.read_csv(f, comment='#', delimiter='\t',header=None,usecols=[0,1,7,9], names=headers)
                #Convert to a dictionary for processing clearly

            rows = reader.T.to_dict().values()

            for row in rows:
                row['Size'], row['Event_Type'],row['End']=parse_event(row)
                # do not include small insertion/deletion calls from Pindel
                if row['Size'] in range(-10,10):
                    continue
                #only normal chr are process, GL is ignored
                try:
                    if chrm:
                        chr1 = 'chr'+str(chromosomes[row['CHROM']])
                    else:
                        chr1 = str(chromosomes[row['CHROM']])
                except KeyError:
                    continue

                #Setup the variables to be returned
                gene1, transcripts1=['Intergenic'], []
                gene2, transcripts2=['Intergenic'], []
                region = 'Intergenic'
                # each segment is assigned to a gene if either the
                # start or end coordinate falls within the feature boundaries.
                start_intervals = exons[chr1].search(int(row['POS']))
                end_intervals = exons[chr1].search(int(row['End']))
                spanning_intervals = exons[chr1].search(int(row['POS']), int(row['End']) + 1)

                if start_intervals:
                    gene1, _, transcripts1 = define_transcripts(start_intervals)
                if end_intervals:
                    gene2, _, transcripts2 = define_transcripts(end_intervals)
                if spanning_intervals:
                    _, regions, _ = define_transcripts(spanning_intervals)
                    if 'EXONIC' in regions:
                        region = 'EXONIC'
                    elif 'UTR' in regions:
                        region = 'UTR'
                    elif 'INTRONIC' in regions:
                        region = 'INTRONIC'

                gene=gene1+gene2
                transcripts=transcripts1+transcripts2

                row['Gene'] =';'.join(str(x) for x in set(gene))
                row['Gene_Region'] = region
                row['Transcripts']=';'.join(str(x) for x in sorted(set(transcripts))) #transcripts #combine_transcripts(set(transcripts)) #
                row['Position']=str(chr1)+':'+str(row['POS'])+'-'+str(row['End'])
                
                if args.multi_reads:
                    row['bbmergedReads']=int(row['bbmergedREADS'].split(',')[-1])
                    row['bwamemReads']=int(row['bwamemREADS'].split(',')[-1])
                else:
                    row['Reads']=int(row['READS'].split(',')[-1])
                
                #show absolute value for size
                row['Size']=abs(row['Size'])
                output.append(row)

    if args.multi_reads:
        output.sort(key=itemgetter('bwamemReads'), reverse=True)  #Sort on reads
        out_fieldnames=['Gene','Gene_Region','Event_Type','Size','Position','bbmergedReads', 'bwamemReads','Transcripts']
    else:
        output.sort(key=itemgetter('Reads'), reverse=True)  #Sort on reads
        out_fieldnames=['Gene','Gene_Region','Event_Type','Size','Position','Reads', 'Transcripts']

    writer = csv.DictWriter(args.outfile, extrasaction='ignore',fieldnames=out_fieldnames, delimiter='\t')
    writer.writeheader()
    writer.writerows(output) 

