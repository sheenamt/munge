import os
import csv
from collections import namedtuple, defaultdict
from operator import itemgetter
import pprint
import logging
import sys
import re
from __init__ import __version__

pfx_pattern = re.compile('(OPX|BRO|MRW|INT|EPI|IMM|IMD|MONC|UNK|TESTDATA)', re.IGNORECASE)
pfx_pattern_old = re.compile('^(OPX|LMG|LMED|CON)', re.IGNORECASE)

log = logging.getLogger(__name__)

def get_location(chr, start, stop, **kwargs):
    """
    Format chr_loc for output (chr:start-end)
    """
    if start == stop:
        chr_loc = 'chr%s:%s' % (chr, start)
    else:
        chr_loc = 'chr%s:%s-%s' % (chr, start, stop)
    return chr_loc


def multi_split(source, splitlist):
    """
    Function to split a string given a string of multiple split points
    """
    output = []
    atsplit = True
    if source is None:
        return None
    else:
        for char in source:
            if char in splitlist:
                atsplit = True
            else:
                if atsplit:
                    output.append(char)
                    atsplit = False
                else:
                    output[-1] = output[-1] + char
    return output


def split_chr_loc(d):

    """
    Function to parse chr_loc(chr:start-end) into chrm, start, end
    """
    output = multi_split(d, 'chr:-')
    #if SNP, there is not end position so set end=start
    try:
        end=output[2].strip()
    except IndexError:
        end=output[1].strip()
    chrm = output[0].strip()
    start = output[1].strip()
    return chrm, start, end


def split_string_in_two(data):
    """
    Return info from one column in two columns
    """
    if not data:
        return '-1', '-1'
    else:
        try:
            output=data.split(',')
            freq=output[0]
            count=output[1]
        except IndexError:
            output=data.split('|')
            freq=output[0]
            count=output[1]
    return freq, count


def build_variant_id(data):
    """
    Construct a variant_id from dict `d` containing keys ....
    """
    d={}
    if data[0]=='Position':
        return 'gendb_link','',''
    else:
        d['chromosome'], d['start'], d['end'] = split_chr_loc(data[0])
        d['ref_base'], d['var_base'] =  data[1], data[2]
        ref_reads=data[15]
        var_reads=data[16]
        variant_id='{chromosome}_{start}_{end}_{ref_base}_{var_base}'.format(**d)
        return variant_id, ref_reads, var_reads

def pfx_ok(pfx, pattern=pfx_pattern):
    """Return True if pfx matches compiled regular expression `pattern`

    """
    return False if pfx is None else bool(pattern.search(pfx))


def fix_pfx(pfx):
    """Normalize pfx (for example, to be used in a database search), but
    only if it looks like a real prefix.

    """
    if pfx_ok(pfx, pattern=pfx_pattern):
        return pfx.replace('-', '_').strip()
    else:
        if pfx_ok(pfx, pattern=pfx_pattern_old):
            return pfx.replace('-', '').strip()
    return pfx.strip()


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


def as_number(val):
    if isinstance(val, int):
        return val
    else:
        return int(''.join(c for c in val if c.isdigit()))

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


Node = namedtuple('Node', 'name start end left right')


def partition(features):
    """Partition `features` into a binary tree where `features` is an
    ordered list of items (name string, start int, end int). The tree
    is constructed of Nodes, each of which is a namedtuple with
    elements 'name start end left right'.

    * name - a string or other object identifying the node
    * start, end - the starting and ending coordinates of a range
    * left - a subtree in which node.start < start for all nodes or None
    * right - a subtree in which node.end > end for all nodes or None

    """

    if features:
        i = len(features) / 2
        left, (name, start, end), right = features[:i], features[i], features[i + 1:]
        assert all(row[1] >= end for row in right)        
        assert all(row[2] <= start for row in left)

        return Node(name, start, end, partition(left), partition(right))


def assign(tree, val):
    """Assign `val` a name if it falls within an interval corresponding to
    a Node in `tree` or None otherwise.

    """

    if tree:
        if val < tree.start:
            return assign(tree.left, val)
        elif val > tree.end:
            return assign(tree.right, val)
        else:
            return tree.name


def check_overlapping(features):
    """Check for elements of `features` with overlapping ranges. In the
    case of overlap, print an informative error message and return
    names and positions of overlapping features.

    """

    features = features[:]
    overlapping = []
    for i in range(len(features)-1):
        prev_name, prev_start, prev_end = features[i]
        name, start, end = features[i+1]
        if prev_end > start:
            overlap = ((prev_name, prev_start, prev_end), (name, start, end))
            overlapping.append(overlap)
            log.warning('overlapping features: ' + pprint.pformat(overlap))
            raise ValueError('overlapping features: ' + pprint.pformat(overlap))
    return overlapping

def get_exons(starts, ends, strand='+'):
    """Return a list of (exon, start, end). `exon` is an int, and `starts`
    and `ends` are strings representing comma-delimited lists of
    chromosomal coordinates. `strand` corresponds to the `strand`
    column of the refGene annotation table; if '-', the order of the
    exon numbers is reversed.

    """
    exon_startpos = [int(i) for i in starts.split(',') if i.strip()]
    exon_endpos = [int(i) for i in ends.split(',') if i.strip()]

    #introns are the sections between exons, but not before or after
    intron_startpos=list([x+1 for x in exon_endpos])
    del intron_startpos[-1]
    intron_endpos=list([x-1 for x in exon_startpos])
    del intron_endpos[0]

    exon_count = range(1, len(exon_startpos) + 1)
    intron_count = range(1, len(intron_startpos) +1)
    
    #reverse the number for naming if on the reverse strand
    if strand == '-':
        exon_count.reverse()
        intron_count.reverse()

    exon_nums=['ex'+str(x) for x in exon_count]
    intron_nums=['int'+str(x) for x in intron_count]

    exons=zip(exon_nums, exon_startpos, exon_endpos)
    introns=zip(intron_nums, intron_startpos, intron_endpos)
    
    return sorted(exons+introns)



def build_trees(file):
    """Returns (genes, exons). Each is a dictionary keyed by chromosome
    (`genes`) or gene (`exons`). Values are trees generated by
    `partition`.

    """
    sort_key = itemgetter(1, 2)  # sorts by (start, end)
    genes = defaultdict(list)
    exon_trees = {}

    for row in read_refgene(file):
        # this fails with a KeyError if a nonstandard chromosome name
        # is encountered (assume these have been filtered out already)

        chr = str(chromosomes[row['chrom']])
        # gene = row['name2']
        # for now, include refgene ID in gene name for validation
        gene = '{name2}:{name}'.format(**row)

        genes[chr].append((gene, int(row['txStart']), int(row['txEnd'])))

        exons = get_exons(row['exonStarts'], row['exonEnds'], row['strand'])
        exons = sorted(exons, key=sort_key)
        check_overlapping(exons)
        exon_trees[gene] = partition(exons)



    # create a gene tree for each chromosome
    gene_trees = {}
    for chr, rows in genes.iteritems():
        rows = sorted(rows, key=sort_key)
        check_overlapping(rows)
        gene_trees[chr] = partition(rows)

    return gene_trees, exon_trees


