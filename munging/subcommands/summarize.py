#!/usr/bin/env python
"""
Description: Given probe reference file and refgene.bed, compute per-gene and summary statistics and output any genes not covered as expected
"""
 
import sys
import subprocess
from csv import DictReader
from os import path
 
def build_parser(parser):
    parser.add_argument('--assay', required=True, help="Assay Reference bed file")
    parser.add_argument('--pref_trans', required=True, help="Gene, RefSeq for assay")
    parser.add_argument('--refgene', required=True, help="Refgene bed file")
    parser.add_argument('--outdir', required=False, help="Output directory for summary scripts")

class exonTracker:
    """
    Keeps track of a gene's exons.  When an interval is inserted, any relevant exons are covered.
    Initially, all exon intervals' values are covered=False
    """
    def __init__(self, exonStarts, exonEnds):
        self.exons = dict(((start,end),False) for start,end in zip(exonStarts,exonEnds))
        assert(len(exonStarts) == len(exonEnds))
    def insert(self, start, end):
        # If a probe begins or ends within an exon, mark that exon as covered.
        for exonStart, exonEnd in self.exons.keys():
            if ((exonStart <= start and end < exonEnd) or
                (exonStart < end and end <= exonEnd)):
                if self.exons[(exonStart,exonEnd)] == False:
                    self.exons[(exonStart,exonEnd)] = True
 
def action(args):
    out = args.outdir if args.outdir else ''
    refseqs = {}
    trans = {}
    refgene_header = ['chrom','chromStart','chromEnd','name', 'refseq','score','strand','thickStart',
                      'thickEnd','itemRgb','exonCount','exonSizes','exonStarts'] 
    probes_header = ['chrom', 'chromStart', 'chromEnd' ]#, 'score', 'strand']
    pref_trans_header = ['Gene', 'RefSeq']

    # 1) Read refGene.bed into the refseqs dictionary
    for line in DictReader(open(args.refgene, 'r'), delimiter='\t', fieldnames=refgene_header):
        refseq = line['refseq']
        name = line['name']
        # Dictionary-ize refgene.bed
        # Insert unseen refseqs into the dictionary; 
        # We asume that refgene only has ONE line per refseq
        if refseq not in refseqs:
            exonStarts = [int(line['chromStart']) + int(exonStart) for exonStart in
                            line['exonStarts'].split(',')]

            exonEnds = [exonStart + int(exonSize) for exonStart, exonSize in
                            zip(exonStarts, line['exonSizes'].split(','))]

            refseqs[refseq] = dict( [('name', name),
                                     ('refseq', line['refseq']),
                                     ('chrom', line['chrom'].strip('chr')),
                                     ('chromStart', int(line['chromStart'])),
                                     ('chromEnd', int(line['chromEnd'])),
                                     ('exonTracker', exonTracker(exonStarts, exonEnds)),
                                     ('bases_covered', 0)])
            
            # Sanity checks
            assert(len(exonStarts) == len(exonEnds))
            for start,end in zip(exonStarts, exonEnds):
                assert(start < end)
                assert(int(line['chromStart']) <= start and start < int(line['chromEnd']))
                assert(int(line['chromStart']) < end and end <= int(line['chromEnd']))
        else:
            sys.stderr.write("Refseq {} is listed twice in refGene!".format(line['refseq']))

    # 2) Using bedtools, calculate how many bases are actually covered for each gene
    # First merge our reference file so each base is only represented once
    # Next, intersect it with refgene to see which bases belong to a gene
    # Finally, also output regions that are not in genes 
    
    merged_probes=path.join(out,'merged_probes.bed')
    write_probes=open(merged_probes, 'w')
    merge_probes_args = ['bedtools', 'merge', '-i', args.assay]
    merge_probes = subprocess.Popen(merge_probes_args, stdout=write_probes) 
    write_probes.close()
    
    intersect_args = ['bedtools', 'intersect', '-wo' ,'-a', merged_probes, '-b', args.refgene]
    intersect = subprocess.Popen(intersect_args, stdout=subprocess.PIPE)

    
    # Parse that output, collecting the number of covered bases per-gene, and annotate refseqs dictionary
    # Note: Communicate returns (stdoutdata, stderrdata), stdout is a giant string, not an iterable
    # Also, the last line is just a newline, which must be skipped
    for line in intersect.communicate()[0].split('\n')[:-1]:
        ls = line.split('\t')
        print(ls)
        refseq = ls[7]          # We pick out the refseq of the gene from refGene that was matched
        overlap = int(ls[-1]) # The '-wo' switch from intersect_args put the amount of overlap here
        print(overlap)
        refseqs[refseq]['bases_covered'] += overlap
        print( refseqs[refseq]['refseq'])
        refseqs[refseq]['exonTracker'].insert(int(ls[1]), int(ls[2]))

        assert(refseqs[refseq]['bases_covered'] <= int(refseqs[refseq]['chromEnd']) - int(refseqs[refseq]['chromStart']))

    # 4) Print per-gene summary
    per_refseq_header = ['gene','refseq','total_bases_targeted','length_of_gene',
                       'fraction_of_gene_covered',
                       'exons_with_any_coverage','total_exons_in_gene']
    per_refseq = open(path.join(out, "per_refseq_summary.txt"), 'w')
    per_refseq.write('\t'.join(per_refseq_header) + '\n')

    # While we're looping through refseqs, count the total bases, exons, and refseqs covered
    total_coding_bases = 0
    total_exons = 0
    gene_count = 0

    for gene in DictReader(open(args.pref_trans, 'r'), delimiter='\t', fieldnames=pref_trans_header):
        transcript = gene['RefSeq'].split('.')[0]
        if transcript.upper()=='REFSEQ':
            continue
        try:
            gene['bases_covered']=refseqs[transcript]['bases_covered']
            #Only count this as a covered gene if it has coverage
            if gene['bases_covered'] > 0:
                gene_count +=1
        except KeyError:
            gene['bases_covered']='RefSeq not found'

        exons = [exon for exon in refseqs[transcript]['exonTracker'].exons.values()].count(True)
        outfields = [gene['Gene'], 
                     gene['RefSeq'],
                     gene['bases_covered'],
                     refseqs[transcript]['chromEnd'] - refseqs[transcript]['chromStart'],
                     round(float(gene['bases_covered']) /
                           float(refseqs[transcript]['chromEnd'] - refseqs[transcript]['chromStart']),3),
                     exons,
                     len(refseqs[transcript]['exonTracker'].exons)]

        per_refseq.write('\t'.join([str(field) for field in outfields]) + '\n')

        total_coding_bases += gene['bases_covered']
        total_exons += exons

    # Calculate total regions covered 
    def calulate_total_covered(probes):
        '''calculate the total regions covered by using the merged probes file'''
        total_cov=0
        with open(probes, 'rU') as p:
            for line in p:
                chrm,start,stop=line.split('\t')
                line_sum=int(stop)-int(start)
                total_cov += line_sum

        return total_cov

    total_bases = calulate_total_covered(merged_probes)

    non_intersect_args = ['bedtools', 'intersect', '-v' ,'-a', merged_probes, '-b', args.refgene]
    non_intersect = subprocess.Popen(non_intersect_args, stdout=subprocess.PIPE)
    
    # 5) Print overall summary
    overall = open(path.join(out, "overall_summary.txt"),'w')

    # Note: The total bases and exon counts are probably slightly overestimated, since refseqs can
    # overlap and share bases.  The number of overlapping bases and exons, however, are neglible
    # and cumbersome to calculate
    overall.write("{} unique bases were targeted\n".format(total_bases))
    overall.write("{} unique coding bases were targeted\n".format(total_coding_bases))
    overall.write("{} unique refseqs had at least one base targeted\n".format(gene_count))
    overall.write("{} total exons had some coverage\n".format(total_exons))

    data=non_intersect.communicate()[0].split('\n')[:-1]
    if data:
        overall.write("The following probes did not intersect with transcription region of any UCSC gene:\n")
        for line in data:
            overall.write(line + "\n")
   
