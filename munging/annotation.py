import os
import csv
from collections import namedtuple, defaultdict
from operator import itemgetter
import pprint
import logging
import sys
import re
from intervaltree import Interval, IntervalTree
from __init__ import __version__
from urllib import urlopen
from StringIO import StringIO as BytesIO
import zlib
from collections import defaultdict
import itertools

pfx_pattern = re.compile('(OPX|BRO|MRW|INT|EPI|IMM|IMD|MONC|UNK|TESTDATA)', re.IGNORECASE)
pfx_pattern_old = re.compile('^(OPX|LMG|LMED|CON)', re.IGNORECASE)

log = logging.getLogger(__name__)

# Various files and data strctures specify chromosomes as strings
# encoding ints, like ('1', '2', ..., 'X'), sometimes as ints (1, 2,
# ... 'X'), and sometimes with a prefix ('chr1', 'chr2', ...,
# 'chrX'). `chromosomes` maps all three to the numeric representation.
chrnums = range(1, 23) + ['X', 'Y']
chromosomes = {'chr{}'.format(c): c for c in chrnums}
chromosomes.update({str(c): c for c in chrnums})
chromosomes.update({c: c for c in chrnums})

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


class UCSCTable(object):
    '''A container class for the parsing functions, used in GenomeIntervalTree.from_table``.'''
    REF_GENE_FIELDS = ['bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames']
    @staticmethod
    def REF_GENE(line):
        return dict(zip(UCSCTable.REF_GENE_FIELDS, line.split(b'\t')))

class IntervalMakers(object):
    '''A container class for interval-making functions, used in GenomeIntervalTree.from_table and GenomeIntervalTree.from_bed.'''


    @staticmethod
    def TX(d):
        return [Interval(int(d['txStart']), int(d['txEnd']), d)]

    @staticmethod
    def CDS(d):
        return [Interval(int(d['cdsStart']), int(d['cdsEnd']), d)]

    @staticmethod
    def EXONS(d):
        exStarts = d['exonStarts'].split(b',')
        exEnds = d['exonEnds'].split(b',')
        intron_count=int(d['exonCount'])-1
        exon_count=int(d['exonCount'])
        strand = d['strand']
        for i in range(exon_count):
            exon_d = d.copy()
            if strand == '+':
                exon_d['exonNum']=str(i+1)
            elif strand == '-':
                exon_d['exonNum']=str(exon_count-i)
            #Since interval trees are not inclusive of upper limit, add one to the exon end boundary
            yield Interval(int(exStarts[i]), int(exEnds[i])+1, exon_d)

            #Setup the intron info
            if i < intron_count:
            #Since interval trees are not inclsive of upper limit, add one to the intron start boundary and not to the end boundary
                intron_start=int(exEnds[i])+1
                intron_end=int(exStarts[i+1])
                new_d=d.copy()
                if new_d['strand']=='-':
                    new_d['intronNum']=str(intron_count - i)
                elif new_d['strand']=='+':
                    new_d['intronNum']=str(i+1)
                yield Interval(intron_start, intron_end, new_d)

def _fix(interval):
    '''
    Helper function for ``GenomeIntervalTree.from_bed and ``.from_table``.

    Data tables may contain intervals with begin >= end. Such intervals lead to infinite recursions and
    other unpleasant behaviour, so something has to be done about them. We 'fix' them by simply setting end = begin+1.
    '''
    if interval.begin >= interval.end:
        log.info("Interval with reversed coordinates (begin >= end) detected when reading data. Interval was automatically fixed to point interval [begin, begin+1).")
        return Interval(interval.begin, interval.begin+1, interval.data)
    else:
        return interval

class GenomeIntervalTree(defaultdict):
    '''
    The data structure maintains a set of IntervalTrees, one for each chromosome.
    It is essentially a ``defaultdict(IntervalTree)`` with a couple of convenience methods
    for reading various data formats.
    '''
    def __init__(self):
        super(GenomeIntervalTree, self).__init__(IntervalTree)

    def addi(self, chrom, begin, end, data=None):
        self[chrom].addi(begin, end, data)

    def __len__(self):
        return sum([len(tree) for tree in self.values()])

    @staticmethod
    def from_table(fileobj=None, url='http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz',
                    parser=UCSCTable.REF_GENE, mode='tx', decompress=None):
        '''
        Index the rows of UCSC tables into a ``GenomeIntervalTree`` 

        The table can be either specified as a ``fileobj`` (in which case the data is read line by line),
        or via an ``url`` (the ``url`` may be to a ``txt`` or ``txt.gz`` file either online or locally).
        The type of the table is specified using the ``parser`` parameter. This is a function that takes a line
        of the file (with no line ending) and returns a dictionary, mapping field names to values. This dictionary will be assigned
        to the ``data`` field of each interval in the resulting tree.

        Finally, there are different ways genes can be mapped into intervals for the sake of indexing as an interval tree.
        One way is to represent each gene via its transcribed region (``txStart``..``txEnd``). Another is to represent using
        coding region (``cdsStart``..``cdsEnd``). Finally, the third possibility is to map each gene into several intervals,
        corresponding to its exons (``exonStarts``..``exonEnds``).

        The mode, in which genes are mapped to intervals is specified via the ``mode`` parameter. The value can be ``tx``, ``cds`` and
        ``exons``, corresponding to the three mentioned possibilities.

        The ``parser`` function must ensure that its output contains the field named ``chrom``, and also fields named ``txStart``/``txEnd`` if ``mode=='tx'``,
        fields ``cdsStart``/``cdsEnd`` if ``mode=='cds'``, and fields ``exonCount``/``exonStarts``/``exonEnds`` if ``mode=='exons'``.

        The ``decompress`` parameter specifies whether the provided file is gzip-compressed.
        This only applies to the situation when the url is given (no decompression is made if fileobj is provided in any case).
        If decompress is None, data is decompressed if the url ends with .gz, otherwise decompress = True forces decompression.

        '''
        #Read in data from URL if file not provided
        if fileobj is None:
            data = urlopen(url).read()
            if (decompress is None and url.endswith('.gz')) or decompress:
                data = zlib.decompress(data, 16+zlib.MAX_WBITS)
            fileobj = BytesIO(data)

        interval_lists = defaultdict(list)

        #Setup the interval type
        if mode == 'tx':
            interval_maker = IntervalMakers.TX
        elif mode == 'cds':
            interval_maker = IntervalMakers.CDS
        elif mode == 'exons':
            interval_maker = IntervalMakers.EXONS
        elif getattr(mode, __call__, None) is None:
            raise Exception("Parameter `mode` may only be 'tx', 'cds', 'exons' or a callable")
        else:
            interval_maker = mode

        #Parse the genome data
        for ln in fileobj:
            if not isinstance(ln, bytes):
                ln = ln.encode()
            if ln.startswith('#'):
                continue
            ln = ln.strip()
            d = parser(ln)
            for interval in interval_maker(d):
                interval_lists[d['chrom']].append(_fix(interval))
                
        # Now convert interval lists into trees
        gtree = GenomeIntervalTree()
        for chrom, lst in getattr(interval_lists, 'iteritems', interval_lists.items)():
            gtree[chrom] = IntervalTree(lst)
        return gtree
        
    def __reduce__(self):
        t = defaultdict.__reduce__(self)
        return (t[0], ()) + t[2:]

def ranges(i):
    for a, b in itertools.groupby(enumerate(i), lambda (x, y): y - x):
        b = list(b)
        yield b[0][1], b[-1][1]


class Transcript(object):
    """Helper class to encapsulate transcript data"""
    
    def __init__(self, gene, accession):
        self.gene = gene
        self.accession = accession
        self.exons = []
        self.introns = []

    def __lt__(self, other):
        self_string = self.gene + self.accession
        other_string = other.gene + other.accession
        return self_string < other_string

        
def define_transcripts(chrm_data):
    """Given the interval, set the gene, region and transcripts"""
    
    gene_list, region_list, transcript_list = [], [], []
    
    # create Transcript dictionary from chrm_data
    transcripts = {}
    for start, stop, data in chrm_data:
        gene = data['name2']
        accession = data['name']

        if not transcripts.has_key(accession):
            transcripts[accession] = Transcript(gene, accession)

        t = transcripts[accession]
        
        if data.has_key('exonNum'):
            #print(accession + ': E-' + data['exonNum'])
            t.exons.append(int(data['exonNum']))

        elif data.has_key('intronNum'):
            #print(accession + ': I-' + data['intronNum'])
            t.introns.append(int(data['intronNum']))

    
    # populate the output lists
    for t in transcripts.values():
        refseq='{}:{}'.format(t.gene,t.accession)

        # if t contains only a single exon
        if len(t.introns) == 0:
            transcript_annotation = '{}(exon {})'.format(refseq, str(t.exons[0]))
            region_annotation = 'EXONIC'

        # if t contains only a single intron
        elif len(t.exons) == 0:
            transcript_annotation = '{}(intron {})'.format(refseq, str(t.introns[0]))
            region_annotation = 'INTRONIC'

        else:
            # can I use first index & last index instead of min & max?
            min_e, min_i = min(t.exons), min(t.introns)
            if min_e <= min_i:
                head = 'exon ' + str(min_e)
            else:
                head = 'intron ' + str(min_i)

            max_e, max_i = max(t.exons), max(t.introns)
            if max_e > max_i:
                tail = 'exon ' + str(max_e)
            else:
                tail = 'intron ' + str(max_i)

            transcript_annotation = '{}({} - {})'.format(refseq, head, tail)
            region_annotation = 'EXONIC-INTRONIC'

        gene_list.append(t.gene)
        region_list.append(region_annotation)
        transcript_list.append(transcript_annotation)

    return sorted(set(gene_list)), sorted(set(region_list)), sorted(set(transcript_list))


