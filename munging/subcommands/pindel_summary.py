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

# def define_transcripts(chrm_data):
#     """Given the interval, set the gene, region and transcripts"""
#     gene1, region, transcript,transcripts=[],[],{},[]
#     for start, stop, data in chrm_data: 
#         gene1.append(data['name2'])
#         refseq='{}:{}'.format(data['name2'],data['name'])
#         if 'exonNum' in data.keys():
#             region.append('Exonic')
#             if transcript.has_key(refseq):
#                 if transcript[refseq].has_key('exons'):
#                     transcript[refseq]['exons'].append(int(data['exonNum']))
#                 else:
#                     transcript[refseq].update({'exons':[int(data['exonNum'])]})
#             else:
#                 transcript[refseq]={'exons':[int(data['exonNum'])]}
#             transcript[refseq]['exons'].sort()
#         if 'intronNum' in data.keys():
#             region.append('Intronic')
#             if transcript.has_key(refseq):
#                 if transcript[refseq].has_key('introns'):
#                     transcript[refseq]['introns'].append(int(data['intronNum']))
#                 else:
#                     transcript[refseq].update({'introns':[int(data['intronNum'])]})
#             else:
#                 transcript[refseq]={'introns':[int(data['intronNum'])]}
#             transcript[refseq]['introns'].sort()

#     for refseq in transcript:
#         if transcript[refseq].has_key('exons'):
#             exon_ranges=list(ranges(transcript[refseq]['exons']))
#             for r in exon_ranges:
#                 if r[0]==r[1]:
#                     transcripts.append('{}(exon {})'.format(refseq,r[0]))
#                 else:
#                     transcripts.append('{}(exons {}-{})'.format(refseq,r[0],r[1]))
#         if transcript[refseq].has_key('introns'):
#             intron_ranges=list(ranges(transcript[refseq]['introns']))
#             for i in intron_ranges:
#                 if i[0]==i[1]:
#                     transcripts.append('{}(intron {})'.format(refseq,i[0]))
#                 else:
#                     transcripts.append('{}(introns {}-{})'.format(refseq,i[0],i[1]))
#     return gene1, region, transcripts
    
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
                gene1, region1, transcripts1=['Intergenic'],['Intergenic'],[]
                gene2, region2, transcripts2=['Intergenic'],['Intergenic'],[]
                # each segment is assigned to a gene if either the
                # start or end coordinate falls within the feature boundaries.
                chrm_start=exons[chr1].search(int(row['POS']))
                chrm_stop=exons[chr1].search(int(row['End']))

                if chrm_start:
                    gene1, region1, transcripts1=define_transcripts(chrm_start)
                if chrm_stop:
                    gene2, region2, transcripts2=define_transcripts(chrm_stop)
                gene=gene1+gene2
                region=region1+region2
                transcripts=transcripts1+transcripts2
                row['Gene'] =';'.join(str(x) for x in set(gene))
                row['Gene_Region']=';'.join(str(x) for x in set(region))
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

