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
from munging.annotation import chromosomes,GenomeIntervalTree, UCSCTable, define_transcripts

csv.field_size_limit(10000000)
log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument('refgene', 
                        help='RefGene file, filtered by preferred transcripts file')
    parser.add_argument('pindel_vcfs', action='append', nargs='+',
                        help='Input files which are vcfs from pindel output')
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
    if abs(size)>1 and data['POS'] == info['END']:
        end=int(info['END'])+size
    else:
        end=int(info['END'])
    return size,svtype, end

def define_transcripts(chrm_data):
    """Given the interval, set the gene, region and transcripts"""
    gene1, region, transcripts=[],[],[]
    for start, stop, data in chrm_data: 
        gene1.append(data['name2'])
        if 'exonNum' in data.keys():
            region.append('Exonic')
            transcript='{}:{}(exon {})'.format(data['name2'],data['name'],data['exonNum'])
            transcripts.append(transcript)
        if 'intronNum' in data.keys():
            region.append('Intronic')
            transcript='{}:{}(intron {})'.format(data['name2'],data['name'],data['intronNum'])
            transcripts.append(transcript)
    return gene1, region, sorted(transcripts)

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
            #Skip the header lines, read in only the columns we need because some unnecessary columns can be millions of characters long
            headers=['CHROM','POS','INFO','READS']
            reader = pandas.read_csv(f, comment='#', delimiter='\t',header=None,usecols=[0,1,7,9], names=headers)
            #Convert to a dictionary for processing clearly
            rows = reader.T.to_dict().values()

            for row in rows:
                row['Size'], row['Event_Type'],row['End']=parse_event(row)
                if row['Size'] in range(-10,10):
                    continue
                #only normal chr are process, GL amd MT are ignored
                try:
                    if chrm:
                        chr1 = 'chr'+str(chromosomes[row['CHROM']])
                    else:
                        chr1 = str(chromosomes[row['CHROM']])
                except KeyError:
                    continue

                #Setup the variables to be returned
                gene1, region, transcripts=[],[],[]
                # each segment is assigned to a gene if either the
                # start or end coordinate falls within the feature boundaries.
                chrm_exons=exons[chr1].search(int(row['POS']), int(row['End']))
                chrm_start=exons[chr1].search(int(row['POS']))
                chrm_stop=exons[chr1].search(int(row['End']))

                #Usual case: both start and stop are in a coding region
                if chrm_exons:
                    gene1, region, transcripts=define_transcripts(chrm_exons)
                # #But if start isn't in coding, but stop is, currently an error
                elif chrm_stop and not chrm_start:
                    print 'hit this case'
                    sys.exit()
                    gene1, region, transcripts=define_transcripts(chrom_stop)
                #Otherwise if neither start nor stop are in coding, label everything as intergenic
                else:
                    gene1=['Intergenic',]
                    region=['Intergenic',]
                    transcripts=[]

                row['Gene'] =';'.join(str(x) for x in set(gene1))
                row['Gene_Region']=';'.join(str(x) for x in set(region))
                row['Transcripts']=';'.join(str(x) for x in set(transcripts))

                row['Position']=str(chr1)+':'+str(row['POS'])+'-'+str(row['End'])
                row['Reads']=int(row['READS'].split(',')[-1])

                output.append(row)

    output.sort(key=itemgetter('Gene'))  #Sort on reads
    output.sort(key=itemgetter('Reads'), reverse=True)  #Sort on reads

    out_fieldnames=['Gene','Gene_Region','Event_Type','Size','Position','Reads', 'Transcripts']
    writer = csv.DictWriter(args.outfile, extrasaction='ignore',fieldnames=out_fieldnames, delimiter='\t')
    writer.writeheader()
    writer.writerows(output) 

