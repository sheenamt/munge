import os
import csv
import pandas as pd
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
import itertools

pfx_pattern = re.compile('(OPX|BRO|MRW|INT|EPI|IMM|IMD|MONC|UNK|TESTDATA)', re.IGNORECASE)
pfx_pattern_old = re.compile('^(OPX|LMG|LMED|CON)', re.IGNORECASE)

log = logging.getLogger(__name__)

# Various files and data strctures specify chromosomes as strings
# encoding ints, like ('1', '2', ..., 'X'), sometimes as ints (1, 2,
# ... 'X'), and sometimes with a prefix ('chr1', 'chr2', ...,
# 'chrX'). `chromosomes` maps all three to the numeric representation.
chrnums = [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT']
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

def ranges(i):
    for a, b in itertools.groupby(enumerate(i), lambda (x, y): y - x):
        b = list(b)
        yield b[0][1], b[-1][1]

class UCSCTable(object):
    '''A container class for the parsing functions, used in GenomeIntervalTree.from_table``.'''
    URL = 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz'
    REF_GENE_FIELDS = ['bin',
                      'name',
                      'chrom',
                      'strand',
                      'txStart',
                      'txEnd',
                      'cdsStart',
                      'cdsEnd',
                      'exonCount',
                      'exonStarts',
                      'exonEnds',
                      'score',
                      'name2',
                      'cdsStartStat',
                      'cdsEndStat',
                      'exonFrames']

    def __init__(self, fileobj=None, decompress=None):
        if fileobj is None:
            data = urlopen(URL).read()
            if (decompress is None and URL.endswith('.gz')) or decompress:
                data = zlib.decompress(data, 16+zlib.MAX_WBITS)
            fileobj = BytesIO(data)

        self.data = pd.read_csv(fileobj, sep='\t', header=None, comment='#', names=UCSCTable.REF_GENE_FIELDS)    

class GenomeIntervalTree(defaultdict):
    '''
    The data structure maintains a set of IntervalTrees, one for each chromosome.
    It is essentially a ``defaultdict(IntervalTree)`` with a couple of convenience methods
    for reading various data formats.
    '''
    def __init__(self):
        super(GenomeIntervalTree, self).__init__(IntervalTree)

    def addi(self, data):
        """
        ADD A DOCSTRING
        """
        t = Transcript(data)
        chrom = t.chrom
        begin = t.tx_start
        end = t.tx_end
        self[chrom].addi(begin, end, t)

    def __len__(self):
        return sum([len(tree) for tree in self.values()])

    def __reduce__(self):
        t = defaultdict.__reduce__(self)
        return (t[0], ()) + t[2:]

    @staticmethod
    def from_table(fileobj=None, decompress=None):
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
        gtree = GenomeIntervalTree()
        table = UCSCTable(fileobj=fileobj, decompress=decompress)

        for _, row in table.data.iterrows():
            gtree.addi(row)
                
        return gtree

class Transcript(object):
    """
    ADD CLASS DOCUMENTATION

    ALL COORDINATES ARE 0-BASED
    UCSC downloaded tables use 0-based start and 1-based end for some incomprehensible reason
    see http://genome.ucsc.edu/FAQ/FAQtracks#tracks1
    """

    @staticmethod
    def _interval_to_str(interval):
        """
        ADD A DOCSTRING
        """
        data = interval[2]
        if data[0] == 'UTR':
            return 'UTR'
        elif data[0] == 'EXON':
            return "exon {}".format(str(data[1]).zfill(2))
        elif data[0] == 'INTRON':
            return "intron {}".format(str(data[1]).zfill(2))
        else:
            raise ValueError("An interval must be 'EXON', 'INTRON', or 'UTR'.")

    def __init__(self, data):
        """
        ADD A DOCSTRING
        """
        # populate fields from data
        self.gene = data['name2']
        self.id = data['name']
        self.chrom = data['chrom']
        if 'chr' in self.chrom:
            self.chrom = self.chrom.replace('chr', '')
        self.strand = data['strand']
        if self.strand not in ['+', '-']:
            raise ValueError("A transcript must be on the '+' or '-' strand")
        self.tx_start = int(data['txStart'])
        self.tx_end = int(data['txEnd']) - 1    # convert end to 0-based, inclusive coordinate
        self.cd_start = int(data['cdsStart'])
        self.cd_end = int(data['cdsEnd']) - 1   # convert end to 0-based, inclusive coordinate
        self.exon_count = int(data['exonCount'])
        # these comma separated entries end with a comma
        # so drop the last entry from the str.split() array
        self.exon_starts = [ int(i) for i in data['exonStarts'].split(',')[0:-1] ]
        self.exon_ends = [ int(i) - 1 for i in data['exonEnds'].split(',')[0:-1] ]  # convert end to 0-based, inclusive coordinate
        self.exon_frames = [ int(i) for i in data['exonFrames'].split(',')[0:-1] ]

        # create IntervalTree for exons, introns, and UTRs
        self.tree_with_utrs = IntervalTree()
        # create IntervalTree without UTRs
        self.tree_no_utrs = IntervalTree()

        # add the UTRs
        # define 5' or 3' based on strand
        if self.strand == '+':
            first_utr = 5
            second_utr = 3
        else:
            first_utr = 3
            second_utr = 5

        # handle edge cases where there is no UTR
        if self.tx_start < self.cd_start:
            self.tree_with_utrs[self.tx_start : self.cd_start] = ('UTR', first_utr)
        if self.tx_end > self.cd_end:
            self.tree_with_utrs[self.cd_end + 1 : self.tx_end + 1] = ('UTR', second_utr)

        # add the exons
        for i in range(self.exon_count):
            exon_start = self.exon_starts[i]
            exon_end = self.exon_ends[i]

            # define exon number based on strand
            if self.strand == '+':
                exon_num = i + 1
            else:
                exon_num = self.exon_count - i
            
            # since interval trees are not inclusive of upper limit, add one to the exon end boundary
            self.tree_no_utrs[exon_start : exon_end + 1] = ('EXON', exon_num)

            # if the exon is entirely encompassed by a UTR, don't add it to the UTR tree
            if exon_end < self.cd_start or exon_start > self.cd_end:
                    continue
            # if the exon is split by a UTR, adjust its boundary
            start = max(exon_start, self.cd_start)
            end = min(exon_end, self.cd_end)
            self.tree_with_utrs[start : end + 1] = ('EXON', exon_num)

        # add the introns
        for i in range(self.exon_count - 1):
            intron_start = self.exon_ends[i] + 1
            intron_end = self.exon_starts[i + 1] - 1

            if self.strand == '+':
                intron_num = i + 1
            else:
                intron_num = self.exon_count - 1 - i
                
            self.tree_no_utrs[intron_start : intron_end + 1] = ('INTRON', intron_num)

            # if the intron is entirely encompassed by a UTR, don't add it to the UTR tree
            if intron_end < self.cd_start or intron_start > self.cd_end:
                    continue
            
            self.tree_with_utrs[intron_start : intron_end + 1] = ('INTRON', intron_num)

    def __str__(self):
        """
        ADD A DOCSTRING
        """
        return "{}:{}".format(self.gene, self.id)

    def __lt__(self, other):
        # change to compare by chr:pos?
        """
        ADD A DOCSTRING
        """
        return str(self) < str(other)

    def __len__(self):
        return self.tx_end - self.tx_start + 1

    def get_annotation(self, start, stop=None, report_utr=True):
        # consider idiot proofing for when start == stop
        """
        ADD A DOCSTRING
        """
        if not stop:
            stop = start + 1
        if report_utr:
            intervals = sorted(self.tree_with_utrs[start:stop])
        else:
            intervals = sorted(self.tree_no_utrs[start:stop])

        # if [start,stop] doesn't overlap with the Transcript, return an empty string
        if len(intervals) == 0:
            return ''

        # otherwise, return <gene>:<id><suffix>
        prefix = str(self)

        # if the search space covers a single interval, suffix is of form (exon 02)
        if len(intervals) == 1:
            single = intervals[0]
            suffix = "({})".format(Transcript._interval_to_str(single))
        # otherwise suffix is of form (exon 02 - intron 09)
        else:
            if self.strand == '+':
                first = intervals[0]
                last = intervals[-1]
            else:
                first = intervals[-1]
                last = intervals[0]

            suffix = ":({} - {})".format(Transcript._interval_to_str(first), Transcript._interval_to_str(last))

        return prefix + suffix

    def get_exons(self, start, stop=None, report_utr=True):
        # consider idiot proofing for when start == stop
        """
        ADD A DOCSTRING
        """
        if not stop:
            stop = start + 1
        if report_utr:
            intervals = self.tree_with_utrs[start:stop]
        else:
            intervals = self.tree_no_utrs[start:stop]

        # return a list of every exon number in the query space
        exons = [i[2][1] for i in intervals if i[2][0] == 'EXON']
        return sorted(exons)

    def get_region_types(self, start, stop=None, report_utr=True):
        # consider idiot proofing for when start == stop
        """
        ADD A DOCSTRING
        """
        if not stop:
            stop = start + 1
        if report_utr:
            intervals = self.tree_with_utrs[start:stop]
        else:
            intervals = self.tree_no_utrs[start:stop]

        # return a list of every unique region type in query space
        region_types = set(i[2][0] for i in intervals)
        return sorted(list(region_types))

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

        if data.has_key('intronNum'):
            #print(accession + ': I-' + data['intronNum'])
            t.introns.append(int(data['intronNum']))

        if data.has_key('UTR'):
            #print(accession + ': U-' + data['UTR'])
            t.utrs.append(data['UTR'])
    
    # populate the output lists
    for t in transcripts.values():
        gene_list.append(t.gene)
        region_list.append(t.get_type())
        transcript_list.append(t.get_annotation())

    return sorted(set(gene_list)), sorted(set(region_list)), sorted(set(transcript_list))

