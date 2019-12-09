import pandas as pd
from collections import namedtuple, defaultdict
from operator import itemgetter
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
# 'chrX'). `chromosomes` maps all three to the string representation.
chrnums = [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT']
chromosomes = {'chr{}'.format(c): c for c in chrnums}
chromosomes.update({str(c): c for c in chrnums})
chromosomes.update({i: str(i) for i in range(1, 23)})

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

def _check_start_stop(start, stop):
    if stop is None:
        stop = start + 1
    if stop <= start:
        raise ValueError("start must be < stop")
    return start, stop

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
        # Intervals are not inclusive of the end point, so increment when adding
        self[chrom].addi(begin, end + 1, t)

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

        The ``decompress`` parameter specifies whether the provided file is gzip-compressed.
        This only applies to the situation when the url is given (no decompression is made if fileobj is provided in any case).
        If decompress is None, data is decompressed if the url ends with .gz, otherwise decompress = True forces decompression.

        '''
        gtree = GenomeIntervalTree()
        table = UCSCTable(fileobj=fileobj, decompress=decompress)
        for row in table.data.to_dict(orient='records'):
            gtree.addi(row)

        return gtree

class SubTranscript(object):
    """Superclass for Exon, Intron, and UTR"""

    def __init__(self, number, start, end):
        self.number = number
        self.start = start
        self.end = end

    def __lt__(self, other):
        """
        Returns True if self is transcribed before other
        
        5' UTR < Exon 1 < Intron 1 < Exon 2 < ... < 3' UTR

        NOTE: only designed to compare elements from the same gene / transcript
        """
        # 5' UTRs are always first, 3' UTRs are always last
        if isinstance(self, UTR):
            if self.number == 5:
                return True
            else:   # it's 3'
                return False
        # 5' UTRs are always first, 3' UTRs are always last
        elif isinstance(other, UTR):
            if other.number == 3:
                return True
            else:   # it's 5'
                return False
        elif isinstance(other, SubTranscript):
            # Exon 1 comes before Intron 1
            if self.number == other.number:
                if isinstance(self, Exon):
                    return True
                else:
                    return False
            # Exon 1 comes before Exon 2 
            else:
                return self.number < other.number
        else:
            raise TypeError("Cannot compare a SubTranscript to a {}".format(other.type))


class Exon(SubTranscript):
    """Helper class to encapsulate an Exon"""

    def __init__(self, number, start, end, frame, cd_start=None, cd_end=None):
        """ADD A DOCSTRING"""
        super(Exon, self).__init__(number, start, end)
        self.frame = frame
        self.cd_start = cd_start
        self.cd_end = cd_end

    def __str__(self):
        """Returns a string representation of the Exon in the form 'exon <number>'"""
        return "exon {}".format(str(self.number).zfill(2))

    def is_coding(self, start, stop=None):
        """
        Returns True if any part of the Exon within [start, stop) contains coding DNA, otherwise returns False
        
        If stop is not provided, queries the range [start, stop + 1)
        """
        if stop is None:
            stop = start + 1

        # if cd_start or cd_end were not initialiazed, the entire exon fell within a UTR and is non-coding
        if  (self.cd_start is None) or (self.cd_end is None): 
            return False
        # if no part of [start, stop) falls within [cd_start, cd_end], then return as non-coding
        elif start > self.cd_end or stop <= self.cd_start:
            return False
        # otherwise return as coding
        else:
            return True


class Intron(SubTranscript):
    """
    Helper class to encapsulate an Intron
    """
    def __init__(self, number, start, end):
        super(Intron, self).__init__(number, start, end)

    def is_coding(self, start, stop=None):
        """Returns False, an Intron is never coding"""
        return False

    def __str__(self):
        """Returns a string representation of the Intron in the form 'intron <number>'"""
        return "intron {}".format(str(self.number).zfill(2))


class UTR(SubTranscript):
    """
    Helper class to encapsulate a UTR
    """
    def __init__(self, number, start, end):
        super(UTR, self).__init__(number, start, end)
    
    def is_coding(self, start, stop=None):
        """Returns False, a UTR is never coding"""
        return False

    def __str__(self):
        """Returns a string representation of the UTR in the form 'UTR'"""
        return "UTR"


class Transcript(object):
    """
    ADD CLASS DOCUMENTATION

    ALL COORDINATES ARE 0-BASED
    UCSC downloaded tables use 0-based start and 1-based end for some incomprehensible reason
    see http://genome.ucsc.edu/FAQ/FAQtracks#tracks1
    
    In this class, all coordinates are 0-based and inclusive ( e.g. for +-stranded transcripts,
    tx_start is the first base transcribed, tx_end is the last base transcribed)

    For class methods such as get_annotation(), start is inclusive but stop is not (which follows the logic of
    querying the underlying IntervalTree)
    """

    def __init__(self, data):
        """
        Creates an instance of a Transcript object from data: a row (pandas Series or dict) from a UCSC RefGene table.
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
        # these comma separated entries end with a comma, so drop the last entry from the str.split() array
        self.exon_starts = [ int(x) for x in data['exonStarts'].split(',')[0:-1] ]
        self.exon_ends = [ int(x) - 1 for x in data['exonEnds'].split(',')[0:-1] ]  # convert end to 0-based, inclusive coordinate
        self.exon_frames = [ int(x) for x in data['exonFrames'].split(',')[0:-1] ]

        # create IntervalTree
        self.tree = IntervalTree()

        # add the UTRs
        # define 5' or 3' based on strand
        if self.strand == '+':
            first_utr = 5
            second_utr = 3
        else:
            first_utr = 3
            second_utr = 5

        # exclude the edge cases where there is no UTR
        if self.tx_start < self.cd_start:
            self.tree[self.tx_start : self.cd_start] = UTR(first_utr, self.tx_start, self.cd_start - 1)
        if self.tx_end > self.cd_end:
            self.tree[self.cd_end + 1 : self.tx_end + 1] = UTR(second_utr, self.cd_end + 1, self.tx_end)

        # add the exons
        for i in range(self.exon_count):
            exon_start = self.exon_starts[i]
            exon_end = self.exon_ends[i]
            exon_frame = self.exon_frames[i]

            # define exon number based on strand
            if self.strand == '+':
                exon_num = i + 1
            else:
                exon_num = self.exon_count - i

            # if the exon is split by a UTR, adjust its boundary
            start = max(exon_start, self.cd_start)
            end = min(exon_end, self.cd_end)

            # since interval trees are not inclusive of upper limit, add one to the exon end boundary
            if start > exon_end or end < exon_start:
                self.tree[exon_start : exon_end + 1] = Exon(exon_num, exon_start, exon_end, exon_frame)
            else:
                self.tree[exon_start : exon_end + 1] = Exon(exon_num, exon_start, exon_end, exon_frame, cd_start=start, cd_end=end)

        # add the introns
        for i in range(self.exon_count - 1):
            intron_start = self.exon_ends[i] + 1
            intron_end = self.exon_starts[i + 1] - 1

            # handle edge cases where two exons are adjacent, with no intervening intron (e.g. ZNF274)
            if intron_start == intron_end + 1:
                continue

            if self.strand == '+':
                intron_num = i + 1
            else:
                intron_num = self.exon_count - 1 - i
                
            self.tree[intron_start : intron_end + 1] = Intron(intron_num, intron_start, intron_end)

    def __str__(self):
        """Returns a string representation of the Transcript in the form '<gene>:<id>'"""
        return "{}:{}".format(self.gene, self.id)

    def __lt__(self, other):
        # change to compare by chr:pos?
        """Returns True if the string representation of self comes alphabetically before that of other"""
        return str(self) < str(other)

    def __len__(self):
        """Returns the distance in nucleotides between the start and end of transcription"""
        return self.tx_end - self.tx_start + 1

    def _remove_utr_ons(self, subtranscripts, start, stop):
        """
        Given a list of Subtranscripts, and a range [start, stop), returns the list of Subtranscripts
        where Introns encompassed by a UTR and Exons that contain no coding DNA in [start, stop) have been removed
        
        NOTE: while the subtranscripts argument is not mutated, the returned list contains references to the original Subtranscript objects
        """
        new_list = []
        for s in subtranscripts:
            if isinstance(s, UTR):
                new_list.append(s)
            elif isinstance(s, Exon) and s.is_coding(start, stop):
                new_list.append(s)
            elif isinstance(s, Intron) and s.end >= self.cd_start and s.start <= self.cd_end:
                new_list.append(s)
        return new_list

    def get_annotation(self, start, stop=None, report_utr=True):
        """
        ADD A DOCSTRING
        """
        start, stop = _check_start_stop(start, stop)
        intervals = self.tree[start:stop]
        subtranscripts = [ x[2] for x in intervals ]

        # if reporting UTR, remove introns and exons that are entrirely encompassed by the UTR within [start, stop)
        if report_utr:
            subtranscripts = sorted(self._remove_utr_ons(subtranscripts, start, stop))
        # if not reporting UTR, remove UTRs from subtranscripts
        else:
            subtranscripts = sorted([ x for x in subtranscripts if not isinstance(x, UTR)])

        # if [start,stop] doesn't overlap with the Transcript, return None
        if len(subtranscripts) == 0:
            return None

        # otherwise, return <gene>:<id><suffix>
        prefix = str(self)

        # if the search space covers a single interval, suffix is of form (exon 02)
        if len(subtranscripts) == 1:
            single = subtranscripts[0]
            suffix = "({})".format(str(single))
        # otherwise suffix is of form (exon 02 - intron 09)
        else:
            first = subtranscripts[0]
            last = subtranscripts[-1]
            suffix = ":({} - {})".format(str(first), str(last))

        return prefix + suffix

    def get_exons(self, start=None, stop=None, report_utr=True):
        # consider idiot proofing for when start == stop
        """
        ADD A DOCSTRING
        """
        if not start:
            start = self.tx_start
            stop = self.tx_end

        start, stop = _check_start_stop(start, stop)
        intervals = self.tree[start:stop]
        exons = [ x[2] for x in intervals if isinstance(x[2], Exon) ]
        
        if report_utr:
            exons = [ x for x in exons if x.is_coding(start, stop)]

        return sorted(exons)

    def get_region_types(self, start, stop=None, report_utr=True):
        # consider idiot proofing for when start == stop
        """
        ADD A DOCSTRING
        """
        start, stop = _check_start_stop(start, stop)

        intervals = self.tree[start:stop]
        subtranscripts = [ x[2] for x in intervals ]

        # if reporting UTR, remove introns and exons that are entrirely encompassed by the UTR within [start, stop)
        if report_utr:
            subtranscripts = sorted(self._remove_utr_ons(subtranscripts, start, stop))
        # if not reporting UTR, remove UTRs from subtranscripts
        else:
            subtranscripts = sorted([ x for x in subtranscripts if not isinstance(x, UTR) ])

        # return a list of every unique region type in query space
        region_types = set()
        for s in subtranscripts:
            if isinstance(s, UTR):
                region_types.add('UTR')
            elif isinstance(s, Exon):
                region_types.add('EXONIC')
            elif isinstance(s, Intron):
                region_types.add('INTRONIC')
            else:
                raise TypeError("A SubTranscript must be of type Exon, Intron, or UTR")

        return region_types

def gene_info_from_transcripts(transcripts):
    """ADD A DOCSTRING"""
    genes = [t.gene for t in transcripts]
    if len(genes) == 0:
        return ['Intergenic']
    else:
        return  sorted(set(genes))

def region_info_from_transcripts(transcripts, start, stop=None, report_utr=True):
    """ADD A DOCSTRING"""
    start, stop = _check_start_stop(start, stop)
    region_types = set()
    for t in transcripts: 
        region_types.update(t.get_region_types(start, stop, report_utr))

    if len(region_types) == 0:
        return ['Intergenic']
    else:
        return sorted(region_types)

def transcript_info_from_transcripts(transcripts, start, stop=None, report_utr=True):
    """ADD A DOCSTRING"""
    start, stop = _check_start_stop(start, stop)
    annotations = [t.get_annotation(start, stop, report_utr) for t in transcripts]
    return sorted(set(annotations))

