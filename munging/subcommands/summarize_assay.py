"""
Given probe reference file, list of preferred transcripts and refgene.bed, 
compute per-refseq and summary statistics and output any genes not covered as expected
Requires refgene data in bed format, use refgene_to_bed script to create
"""
 
import sys
import subprocess
from csv import DictReader, DictWriter
import os

from munging.annotation import chromosomes,GenomeIntervalTree, UCSCTable

def build_parser(parser):
    parser.add_argument('--assay', required=True, help="Assay Reference bed file, sorted, with ^M removed from end of lines")
    parser.add_argument('--pref_trans', required=True, help="Gene, RefSeq for assay")
    parser.add_argument('--refgene', required=True, help="UCSC Refgene data in bed format ")
    parser.add_argument('--bedtools', default='',help='Location of bedtools singularity image')
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
        #Probes are merged, and therefore may cover multiple exons
        #probe: 1-100, exons: 1-20, 30-49, 59-99
        for exonStart, exonEnd in self.exons.keys():
            if (int(start) >= int(exonStart) and int(end) < int(exonEnd)) or (int(end) > int(exonStart) and int(end) <= int(exonEnd)) or (int(exonStart) > int(start) and int(exonEnd) < int(end)):
                self.exons[(exonStart,exonEnd)] = True
                
def action(args):
    BEDTOOLS='singularity exec --bind {} --pwd {} {}'.format(os.getcwd(), os.getcwd(), args.bedtools)

    out = args.outdir if args.outdir else ''
    refseqs = {}
    pref_trans = {}

    refgene_header = ['chrom','chromStart','chromEnd','name', 'refseq','exonCount','exonSizes','exonStarts','exonEnds'] 
    probes_header = ['chrom', 'chromStart', 'chromEnd' ]
    pref_trans_header = ['Gene', 'RefSeq']

    
    # 1) Read refGene.txt into the refseqs dictionary
    for line in DictReader(open(args.refgene, 'r'), delimiter='\t', fieldnames=refgene_header):
        refseq = line['refseq']
        name = line['name']
        # Dictionary-ize refgene.bed
        # Insert unseen refseqs into the dictionary; 
        # We asume that refgene only has ONE line per refseq
        exonStarts=[x for x in line['exonStarts'].split(',')]
        exonEnds=[x for x in line['exonEnds'].split(',')]
        if refseq not in refseqs:
            refseqs[refseq] = dict( [('name', name),
                                     ('refseq', line['refseq']),
                                     ('chrom', line['chrom'].strip('chr')),
                                     ('chromStart', int(line['chromStart'])),
                                     ('chromEnd', int(line['chromEnd'])),
                                     ('exonTracker', exonTracker(exonStarts, exonEnds)),
                                     ('bases_covered', 0)])
            
            #Sanity checks
            assert(len(exonStarts) == len(exonEnds))
            for start,end in zip(exonStarts, exonEnds):
                assert(int(start) < int(end))
                assert(int(line['chromStart']) <= int(start) and int(start) < int(line['chromEnd']))
                assert(int(line['chromStart']) < int(end) and int(end) <= int(line['chromEnd']))
        else:
            sys.stderr.write("Refseq {} is listed twice in refGene!".format(line['refseq']))

    # 2) Using bedtools, calculate how many bases are actually covered for each gene
    # First merge our reference file so each base is only represented once
    # Next, intersect it with refgene to see which bases belong to a gene
    # Finally, also output regions that are not in genes 
    
    merged_probes=os.path.join(out,'merged_probes.bed')
    write_probes=open(merged_probes, 'w')
    merge_probes_args = [x for x in BEDTOOLS.split(' ')]+['bedtools', 'merge', '-i', args.assay]
    merge_probes = subprocess.call(merge_probes_args, stdout=write_probes) 
    write_probes.close()

    intersection=os.path.join(out,'intersect_probes_refgene.txt')
    write_intersect=open(intersection, 'w')
    intersect_args = [x for x in BEDTOOLS.split(' ')]+['bedtools','intersect', '-wo' ,'-a', merged_probes, '-b', args.refgene]
    intersect = subprocess.call(intersect_args, stdout=write_intersect)
    write_intersect.close()

    # Parse that output, collecting the number of covered bases per-gene, and annotate refseqs dictionary
    # Note: Communicate returns (stdoutdata, stderrdata), stdout is a giant string, not an iterable
    # Also, the last line is just a newline, which must be skipped
    for line in open(intersection):
        ls = line.split('\t')
        refseq = ls[7]          # We pick out the refseq of the gene from refGene that was matched
        overlap = int(ls[-1]) # The '-wo' switch from intersect_args put the amount of overlap here
        refseqs[refseq]['bases_covered'] += overlap
        refseqs[refseq]['exonTracker'].insert(int(ls[1]), int(ls[2]))

    # 4) Print per-refseq summary
    per_refseq_header = ['gene','refseq','total_bases_targeted','length_of_gene','fraction_of_gene_covered'] 
#                         'fraction_of_gene_covered','exons_with_any_coverage','total_exons_in_gene']
    per_refseq_writer = DictWriter(open(os.path.join(out, "per_refseq_summary.txt"), 'w'), fieldnames=per_refseq_header,  delimiter='\t', extrasaction='ignore')
    per_refseq_writer.writeheader()
    # While we're looping through refseqs, count the total bases, exons, and refseqs covered
    total_coding_bases = 0
    total_exons = 0
    gene_count = 0

    for gene in DictReader(open(args.pref_trans, 'r'), delimiter='\t', fieldnames=pref_trans_header):
        transcript = gene['RefSeq'].split('.')[0]
        if transcript.upper()=='REFSEQ':
            continue
        try:
            #assert that the NM provided in the Pref_Trans is for the Gene they requested
            if refseqs[transcript]['name']!=gene['Gene']:
                outfields = dict([('gene', gene['Gene']), 
                                  ('refseq', gene['RefSeq']),
                                  ('total_bases_targeted', 'Incorrect RefSeq for this Gene'),
                                  ('length_of_gene','NA'),
                                  ('fraction_of_gene_covered','NA'),
                                  ('exons_with_any_coverage','NA'),
                                  ('total_exons_in_gene','NA')])
            else:
                gene['bases_covered']=refseqs[transcript]['bases_covered']
                #Only count this as a covered gene if it has coverage
                if gene['bases_covered'] > 0:
                    gene_count +=1
                exons = [exon for exon in refseqs[transcript]['exonTracker'].exons.values()].count(True)

                outfields = dict([('gene', gene['Gene']), 
                                  ('refseq', gene['RefSeq']),
                                  ('total_bases_targeted', gene['bases_covered']),
                                  ('length_of_gene',refseqs[transcript]['chromEnd'] - refseqs[transcript]['chromStart']),
                                  ('fraction_of_gene_covered',round(float(gene['bases_covered']) /
                                                                    float(refseqs[transcript]['chromEnd'] - refseqs[transcript]['chromStart']),3)),
                                  ('exons_with_any_coverage',exons),
                                  ('total_exons_in_gene',len(refseqs[transcript]['exonTracker'].exons))])
                total_coding_bases += gene['bases_covered']
                total_exons += exons

        #If this refseq isn't found, we should state that, cleanly 
        except KeyError:
            outfields = dict([('gene', gene['Gene']), 
                              ('refseq', gene['RefSeq']),
                              ('total_bases_targeted', 'RefSeq not found'),
                              ('length_of_gene','NA'),
                              ('fraction_of_gene_covered','NA'),
                              ('exons_with_any_coverage','NA'),
                              ('total_exons_in_gene','NA')])
        pref_trans[gene['Gene']] = outfields

    
    #Process all data in the refseq file that overlaps our probes
    for transcript,data in refseqs.iteritems():
        #add data for probes we didn't plan to cover
        if data['name'] in pref_trans.keys() :
            continue
        else:
            if data['bases_covered'] > 0:
                gene_count +=1
                exons = [exon for exon in data['exonTracker'].exons.values()].count(True)
                outfields = dict([('gene', data['name']), 
                                  ('refseq', transcript),
                                  ('total_bases_targeted', data['bases_covered']),
                                  ('length_of_gene',data['chromEnd'] - data['chromStart']),
                                  ('fraction_of_gene_covered',round(float(data['bases_covered']) /
                                                                    float(data['chromEnd'] - data['chromStart']),3)),
                                  ('exons_with_any_coverage',exons),
                                  ('total_exons_in_gene',len(data['exonTracker'].exons))])
                total_coding_bases += data['bases_covered']
                total_exons += exons
                pref_trans[data['name']] = outfields

    for gene,data in sorted(pref_trans.iteritems()):
        per_refseq_writer.writerow(data)



    #5)  Calculate total regions covered 
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
    non_intersection=os.path.join(out,'non_intersect_probes_refgene.txt')
    write_non_intersect=open(non_intersection, 'w')
    non_intersect_args = [x for x in BEDTOOLS.split(' ')]+['bedtools','intersect', '-v' ,'-a', merged_probes, '-b', args.refgene]
    non_intersect = subprocess.call(non_intersect_args, stdout=write_non_intersect)
    write_non_intersect.close()

    # 6) Print overall summary
    overall = open(os.path.join(out, "overall_summary.txt"),'w')

    # Note: The total bases and exon counts are probably slightly overestimated, since refseqs can
    # overlap and share bases.  The number of overlapping bases and exons, however, are neglible
    # and cumbersome to calculate
    overall.write("{} unique bases were targeted\n".format(total_bases))
    overall.write("{} unique bases within gene boundaries were targeted\n".format(total_coding_bases))
    overall.write("{} unique refseqs had at least one base targeted\n".format(gene_count))
#    overall.write("{} total exons had some coverage\n".format(total_exons))

    data=open(non_intersection)
    if data:
        overall.write("The following probes did not intersect with transcription region of any UCSC gene:\n")
        for line in data:
            overall.write(line)
   
