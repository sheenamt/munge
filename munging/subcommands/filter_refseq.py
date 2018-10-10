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

import os
import sys
import csv
from itertools import chain, groupby
from operator import itemgetter
import logging

#from munging import package_data
from munging.utils import Opener
from munging.annotation import (read_refgene, chromosomes, check_overlapping,
                               get_preferred_transcripts)

log = logging.getLogger(__name__)


def get_genes_files_assay_names():
    genes_files = package_data(None, 'genes/*.genes')
    return sorted({os.path.basename(fn).split('_')[0] for fn in genes_files})


def build_parser(parser):
    parser.add_argument('annotations', type=Opener(), metavar='FILE',
                        help='Annotations file')
    parser.add_argument('genes', type=Opener(), metavar='FILE',
                        nargs='+', help="""One or more files defining
                        preferred transcripts""")
    parser.add_argument('-o', '--outfile', type=Opener('w+'), metavar='FILE',
                        default=sys.stdout, help='output file')


def action(args):
    transcripts = chain.from_iterable(csv.reader(file, delimiter='\t')
                                      for file in args.genes)

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
