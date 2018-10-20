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

from munging.utils import Opener
from munging.annotation import chromosomes,GenomeIntervalTree, UCSCTable


log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument('refgene', 
                        help='RefGene file, filtered by preferred transcripts file')
    parser.add_argument('pindel_vcfs', action='append', nargs='+',
                        help='Input files which are vcfs from pindel output')
    parser.add_argument('-o', '--outfile', type=Opener('w'), metavar='FILE',
                        default=sys.stdout, help='output file')


def parse_event(data):
    ''' Return the length and type of event '''
    
    #Parse read depth and SVtype
    info=dict(item.split('=') for item in data['INFO'].split(";") if "=" in item)
    size=int(info['SVLEN'])
    if info['SVTYPE'] == 'RPL':
        svtype='DEL'
    else:
        svtype=info['SVTYPE']
    #Pindel reports insertions wrong by not setting the end position correctly. It may do it with other data, so test based on reported size rather than svtype
    if size >1 and data['POS'] == info['END']:
        end=int(info['END'])+size
    else:
        end=info['END']
    return size,svtype, end

def action(args):
    exons = GenomeIntervalTree.from_table(open(args.refgene, 'r'), parser=UCSCTable.REF_GENE, mode='exons')
    output = []

    #Skip the header lines 
    (pindel_vcfs,) = args.pindel_vcfs
    for vcf in pindel_vcfs:
        with open(vcf, 'rU')  as f:
            # read in the entire input file so that we can sort it
            fieldnames=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','READS']
            reader = csv.DictReader(filter(lambda row: row[0]!='#', f), delimiter='\t',fieldnames=fieldnames)

            rows = list(reader)
            for row in rows:
                row['Size'], row['Event_Type'],row['End']=parse_event(row)
                if row['Size'] in range(-10,10):
                    continue
                # each segment is assigned to a gene if either the
                # start or end coordinate falls within the feature boundaries.
                try:
                    chr1 = 'chr'+str(chromosomes[row['CHROM']])
                except KeyError:
                    print('chrm not being processed: {}'.format(row['CHROM']))
                    continue

                #Setup the variable to be returned
                gene1=[]
                region=''
                transcripts=[]
                
                #Since ranges are inclusive of the lower limit, but non-inclusive of the upper limit,
                #Make sure we cover everything
                #Usual case: both start and stop are in a coding region
                chrm_start=exons[chr1].search(int(row['POS']))
                chrm_stop=exons[chr1].search(int(row['End']))
                chrm_exons=exons[chr1].search(int(row['POS']), int(row['End']))
                if chrm_exons:
                    for start, stop, data in chrm_exons:
                        gene1.append(data['name2'])
                        if 'exonNum' in data.keys():
                            region='Exonic'
                            transcript='{}:{}(exon {})'.format(data['name2'],data['name'],data['exonNum'])
                            transcripts.append(transcript)
                        if 'intronNum' in data.keys():
                            region='Intronic'
                            transcript='{}:{}(intron {})'.format(data['name2'],data['name'],data['intronNum'])
                            transcripts.append(transcript)

                #But if start isn't in coding, but stop is, process stop
                elif not chrm_start.issubset(chrm_exons):
                    if chrm_stop.issubset(chrm_exons):
                        for start, stop, data in chrm_stop:
                            gene1.append(data['name2'])
                            if 'exonNum' in data.keys():
                                region='Exonic'
                                transcript='{}:{}(exon {})'.format(data['name2'],data['name'],data['exonNum'])
                                transcripts.append(transcript)
                            if 'intronNum' in data.keys():
                                region='Intronic'
                                transcript='{}:{}(intron {})'.format(data['name2'],data['name'],data['intronNum'])
                                transcripts.append(transcript)
                #Otherwise if neither start nor stop are in coding, label everything as intergenic
                else:
                    gene1=['Intergenic',]
                    region='Intergenic'
                    transcripts=[]

                row['Gene'] =';'.join(str(x) for x in set(gene1))
                row['Gene_Region']=region
                row['Transcripts']=';'.join(str(x) for x in set(transcripts))

                row['Position']=str(chr1)+':'+str(row['POS'])+'-'+str(row['End'])
                row['Reads']=int(row['READS'].split(',')[-1])

                output.append(row)

    sorted_output = sorted(output, key=itemgetter('Reads'), reverse=True)  #Sort on reads

    out_fieldnames=['Gene','Gene_Region','Event_Type','Size','Position','Reads', 'Transcripts']
    writer = csv.DictWriter(args.outfile, extrasaction='ignore',fieldnames=out_fieldnames, delimiter='\t')
    writer.writeheader()
    writer.writerows(sorted_output) 

