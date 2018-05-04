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

#    merge_probes_args = ['bedtools', 'merge', '-i', args.assay, '-c', '4,5', '-o', 'distinct,distinct' ]
    merge_probes_args = ['bedtools', 'merge', '-i', args.assay]
    merge_probes = subprocess.Popen(merge_probes_args, stdout=subprocess.PIPE) 
    intersect_args = ['bedtools', 'intersect', '-wo', '-a', 'stdin', '-b', args.refgene]
    intersect = subprocess.Popen(intersect_args, stdin=merge_probes.stdout, stdout=subprocess.PIPE)

    # Parse that output, collecting the number of covered bases per-gene, and annotate refseqs dictionary
    # Note: Communicate returns (stdoutdata, stderrdata), stdout is a giant string, not an iterable
    # Also, the last line is just a newline, which must be skipped
    for line in intersect.communicate()[0].split('\n')[:-1]:
        ls = line.split('\t')
        refseq = ls[7]          # We pick out the refseq of the gene from refGene that was matched
        overlap = int(ls[-1]) # The '-wo' switch from intersect_args put the amount of overlap here
        refseqs[refseq]['bases_covered'] += overlap
        refseqs[refseq]['exonTracker'].insert(int(ls[1]), int(ls[2]))

        assert(refseqs[refseq]['bases_covered'] <= int(refseqs[refseq]['chromEnd']) - int(refseqs[refseq]['chromStart']))

    # 3) Now we count exons
    # refseqs now has info on each probe, arranged by refseq number
    # for key, value in refseqs.items():
    #     if value['bases_covered'] > 0:
    #         print(key, value)


    # # Insert each probe interval into each gene.  The exonTracker knows which exons are and aren't
    # # yet covered.
    # for line in DictReader(open(args.assay, 'r'), delimiter='\t', fieldnames=probes_header):
    #     for name in line['name'].split(','):
    #         if name not in refseqs.keys():
    #             print 'gene not found:', name
    #         else:
    #             refseqs[name]['exonTracker'].insert(int(line['chromStart']), int(line['chromEnd']))

    # 4) Print per-gene summary
    per_refseq_header = ['gene','refseq','total_bases_targeted','length_of_gene',
                       'fraction_of_gene_covered',
                       'exons_with_any_coverage','total_exons_in_gene']
    per_refseq = open(path.join(out, "per_refseq_summary.csv"), 'w')
    per_refseq.write('\t'.join(per_refseq_header) + '\n')

    # While we're looping through refseqs, count the total bases, exons, and refseqs covered
    total_bases = 0
    total_exons = 0
    gene_count = 0

    for gene in refseqs.values():
        # Skip refseqs that weren't targeted
        if gene['bases_covered'] == 0:
            continue
        exons = [exon for exon in gene['exonTracker'].exons.values()].count(True)
        outfields = [gene['name'], 
                     gene['refseq'],
                     gene['bases_covered'],
                     gene['chromEnd'] - gene['chromStart'],
                     round(float(gene['bases_covered']) /
                           float(gene['chromEnd'] - gene['chromStart']),2),
                     exons,
                     len(gene['exonTracker'].exons)]

        per_refseq.write('\t'.join([str(field) for field in outfields]) + '\n')
        gene_count += 1
        total_bases += gene['bases_covered']
        total_exons += exons

    # 5) Print overall summary
    overall = open(path.join(out, "overall_summary.csv"),'w')

    # Note: The total bases and exon counts are probably slightly overestimated, since refseqs can
    # overlap and share bases.  The number of overlapping bases and exons, however, are neglible
    # and cumbersome to calculate
    overall.write("{} unique bases were targeted\n".format(total_bases))
    overall.write("{} unique refseq refseqs had at least one base targeted\n".format(gene_count))
    overall.write("{} total exons had some coverage\n".format(total_exons))

    # 6) Print coverage of expected genes
    coverage_header = ['Gene', 'RefSeq', 'Coverage']
    coverage = open(path.join(out, "expected_coverage.csv"),'w')

    coverage.write('\t'.join(coverage_header) + '\n')
    for gene in DictReader(open(args.pref_trans, 'r'), delimiter='\t', fieldnames=pref_trans_header):
        transcript = gene['RefSeq'].split('.')[0]
        if transcript.upper()=='REFSEQ':
            continue
        try:
            gene['bases_covered']=refseqs[transcript]['bases_covered']
            
        except KeyError:
            gene['bases_covered']='RefSeq not found'

        outfields = [gene['Gene'], 
                     gene['RefSeq'],
                     gene['bases_covered']]
        coverage.write('\t'.join([str(field) for field in outfields]) + '\n')


        # if int(refseqs[transcript]['bases_covered'])=0:
        #     print("not covered", line)
        
