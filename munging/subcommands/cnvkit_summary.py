"""Annotate CNVkit cnr output with genes, exons, and other features.
"""
from natsort import natsorted

import sys
import argparse
import csv
from collections import defaultdict
from operator import itemgetter
import logging
import pandas
from munging.utils import Opener
from munging.annotation import chromosomes,GenomeIntervalTree, UCSCTable

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('refgene', 
                        help='RefGene file, filtered by preferred transcripts file')
    parser.add_argument('cnvkit_cnr', type=Opener(),
                        help='*cnr file from CNVkit')
    parser.add_argument('-p', '--plot', type=Opener('w'), metavar='FILE',
                        default=sys.stdout, help='Create output file required for R plotting scripts')
    parser.add_argument('-o', '--outfile', type=Opener('w'), metavar='FILE',
                        default=sys.stdout, help='output file')


def define_transcripts(chrm_data):
    """Given the interval, set the transcripts"""
    gene,transcripts=[],[]
    exon_number=''
    for start, stop, data in chrm_data: 
        gene='{}:{}'.format(data['name2'], data['name'])
        if 'exonNum' in data.keys():
            transcript='{}:(exon {})'.format(gene,data['exonNum'])
            transcripts.append(transcript)
            exon_number=data['exonNum']
        if 'intronNum' in data.keys():
            transcript='{}:(intron {})'.format(gene,data['intronNum'])
            transcripts.append(transcript)
    trans = ';'.join(str(x) for x in natsorted(set(transcripts)))
    return gene, trans, exon_number

def action(args):
    #Create interval tree of introns and exons,  grouped by chr
    exons = GenomeIntervalTree.from_table(open(args.refgene, 'r'), parser=UCSCTable.REF_GENE, mode='exons')
    output = []

    #try to match chr string in reference file and input data
    chrm=False
    if 'chr' in exons.keys()[0]:
        chrm = True

    #Skip the header lines, read in only the columns we need because some unnecessary columns can be millions of characters long
    headers=['Chr','OriStCoordinate','end','Gene','Adjusted.Mean.of.LogRatio','Depth','Weight']
    reader = csv.DictReader(args.cnvkit_cnr, delimiter='\t')
    #Reset the headers to the names we want (changing capitalization)
    reader.fieldnames=headers

    rows = [x for x in reader]

    for row in rows:
        #Only process on target
        if row['Gene']=='Antitarget' and float(row['Depth']) < 10:
            continue 

        #only normal chr are process, GL amd MT are ignored
        try:
            if chrm:
                chr1 = 'chr'+str(chromosomes[row['Chr']])
            else:
                chr1 = str(chromosomes[row['Chr']])
        except KeyError:
            continue

        # each segment is assigned to a gene if either the
        # start or end coordinate falls within the feature boundaries.
        chrm_exons=exons[chr1].search(int(row['OriStCoordinate']), int(row['end']))
        chrm_start=exons[chr1].search(int(row['OriStCoordinate']))
        chrm_stop=exons[chr1].search(int(row['end']))
                
        #Usual case: both start and stop are in a coding region
        if chrm_exons:
            gene,row['Transcripts'],row['Exon.Number']=define_transcripts(chrm_exons)
        # #But if start isn't in coding, but stop is, currently an error
        elif chrm_stop and not chrm_start:
            print 'hit this case'
            sys.exit()
            gene,row['Transcripts'],row['Exon.Number']=define_transcripts(chrom_stop)
            #Otherwise if neither start nor stop are in coding, label everything as intergenic
        else:
            continue

        row['Gene.Sym']=gene

        row['Position']=str(chr1)+':'+str(row['OriStCoordinate'])+'-'+str(row['end'])
        output.append(row)

    sorted_output = natsorted(output, key=itemgetter('Position'))  #Sort by gene name
    out_fieldnames=['Gene','Position','Depth','Adjusted.Mean.of.LogRatio','Weight', 'Transcripts']
    writer = csv.DictWriter(args.outfile, extrasaction='ignore',fieldnames=out_fieldnames, delimiter='\t')
    writer.writeheader()
    writer.writerows(sorted_output) 

    if args.plot:
        out_fieldnames=['Adjusted.Mean.of.LogRatio','Exon.Number','Chr','OriStCoordinate', 'Gene.Sym']
        writer = csv.DictWriter(args.plot, extrasaction='ignore',fieldnames=out_fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(sorted_output) 


