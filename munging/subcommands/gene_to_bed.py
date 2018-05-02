#!/usr/bin/env python
"""
File: gene_to_bed.py
Modified by: Sheena Todhunter

Description: Given probe reference file and refgene.bed, compute per-gene and summary statistics
"""
import os
import csv
from itertools import chain, groupby
from operator import itemgetter
import logging
 
import sys
import subprocess
import argparse

log = logging.getLogger(__name__)
 
def build_parser(parser):
    parser = argparse.ArgumentParser(description="""Convert UCSC refgene.txt to BED
    format""", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--genes', required=True, help="List of genes and preferred transcripts")
    parser.add_argument('--refgene', required=True, help="Refgene bed file")
    parser.add_argument('--outdir', required=False, help="Output directory for summary scripts")

    return parser.parse_args()

refgene_fields = """
bin
name
chrom
strand
txStart
txEnd
cdsStart
cdsEnd
exonCount
exonStarts
exonEnds
score
name2
cdsStartStat
cdsEndStat
exonFrames
""".split()

# Various files and data strctures specify chromosomes as strings
# encoding ints, like ('1', '2', ..., 'X'), sometimes as ints (1, 2,
# ... 'X'), and sometimes with a prefix ('chr1', 'chr2', ...,
# 'chrX'). `chromosomes` maps all three to the numeric representation.
chrnums = range(1, 23) + ['X', 'Y']
chromosomes = {'chr{}'.format(c): c for c in chrnums}
chromosomes.update({str(c): c for c in chrnums})
chromosomes.update({c: c for c in chrnums})


def read_refgene(file):
    """Read open file-like object `file` containing annotations and
    return a sequence of dicts.

    The `annotations` file is available from the UCSC Genome Browser
    website:
    http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz

    To view the schema, go to https://genome.ucsc.edu/cgi-bin/hgTables
    --> choose track "refSeq Genes" --> choose table "refGene" -->
    "describe table schema"

    """

    return csv.DictReader(file, fieldnames=refgene_fields, delimiter='\t')


def get_preferred_transcripts(transcripts):
    """Return a dict of {gene: {set of transcripts}} given `rows`, a
    sequence of tuples in the format ("gene",
    "transcript1[/transcript/...]")

    """

    tdict = defaultdict(set)
    for row in transcripts:
        # file may have more than one column, or only one: skip genes
        # without a preferred transcript.
        try:
            gene, transcript = row[0], row[1]
        except IndexError:
            continue

        # "transcript" may identify more than one transcript separated
        # with a slash
        tdict[gene].update({t.strip().split('.')[0]
                            for t in transcript.split('/') if t.strip()})

    # convert to a plain dictionary to prevent inadvertently adding keys
    return dict(tdict)


"""Filter a file containing the refGene annotation table, limiting to
preferred transcripts.

The `annotations` file is available from the UCSC Genome Browser website:
http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz

To view the schema, go to https://genome.ucsc.edu/cgi-bin/hgTables -->
choose track "refSeq Genes" --> choose table "refGene" --> "describe
table schema"

`genes` is a tab-delimited file with columns "Gene",
"RefSeq" (and maybe more).

When filtering, choose the preferred transcript. If more than one
preferred transcript is defined, choose one arbitrarily. Failure to
find the preferred transcript for each gene results in an error.

Overlapping genes will result in an error.

"""



# def build_parser(parser):
#     parser.add_argument('annotations', type=Opener(), metavar='FILE',
#                         help='Annotations file')
#     parser.add_argument('genes', type=Opener(), metavar='FILE',
#                         nargs='+', help="""One or more files defining
#                         preferred transcripts""")
#     parser.add_argument('-o', '--outfile', type=Opener('w+'), metavar='FILE',
  #                         default=sys.stdout, help='output file')


def main(args):
    print('ker')
    print(args.genes)
    transcripts = chain.from_iterable(csv.reader(file, delimiter='\t')
                                      for file in args.genes)

    print(transcripts)

    tdict = get_preferred_transcripts(transcripts)

    # read and filter the refgene file
    refseqs = read_refgene(args.annotations)
    fieldnames = refseqs.fieldnames

    refseqs = [r for r in refseqs
               if r['chrom'] in chromosomes and r['name2'] in tdict]

    # sort by chromosome, transcription start
    refseqs.sort(key=lambda row: (
        chromosomes[row['chrom']], int(row['txStart'])))

    # group by gene and choose one transcript for each
    output = []
    for gene, grp in groupby(refseqs, itemgetter('name2')):
        grp = list(grp)
        preferred = tdict[gene]

        if preferred:
            keep = [r for r in grp if r['name'] in preferred]
            if not keep:
                log.error('Error: %s has a preferred transcript of %s but only %s was found' %
                          (gene, preferred, ','.join(r['name'] for r in grp)))
                sys.exit(1)
            elif len(keep) > 1:
                log.warning('{} has more than one preferred transcript; using {}'.format(
                    gene, keep[0]['name']))
        else:
            log.warning('no preferred transcript for {}'.format(gene))
            keep = grp[:1]

        output.append(keep[0])

    # all transcripts are found among preferred transcripts
    assert all(k['name'] in reduce(set.union, tdict.values()) for k in output)

    # Check for overlapping genes and exit with an error if any are
    # found.
    overlapping = []
    rows = sorted(output, key=itemgetter('chrom'))
    for chrom, grp in groupby(rows, key=itemgetter('chrom')):
        grp = sorted(
            [(row['name2'], int(row['txStart']), int(row['txEnd'])) for row in grp],
            key=itemgetter(1, 2)
        )
        overlapping.extend(check_overlapping(grp))

    if overlapping:
        log.error('Error: overlapping genes were found')
        sys.exit(1)

    writer = csv.DictWriter(args.outfile, fieldnames=fieldnames, delimiter='\t')
    writer.writerows(output)

if __name__ == '__main__':
    sys.exit(main())

