"""
Test the pindel_summary script
"""

import subprocess
import filecmp
import logging
import os

from munging.subcommands import pindel_summary
from intervaltree import Interval
from __init__ import TestBase
import __init__ as config

log = logging.getLogger(__name__)


pindel_testfiles = os.path.join(config.datadir, 'pindel')

class TestPindelSummary(TestBase):
    """
    Test the pindel_summary script, which combines vcfs and annotates
    """
    def setUp(self):
        self.outdir = self.mkoutdir()
        self.refgene = os.path.join(pindel_testfiles, 'refgene_test.txt')
        self.data=({'INFO': 'END=1808032;HOMLEN=0;SVLEN=-32;SVTYPE=RPL;NTLEN=32', 'FORMAT': 'GT:AD', 'CHROM': '4', 'POS': '1808000', 
                    'FILTER': 'PASS', 'QUAL': '.', 'READS': '0/0:409,6', 'ALT': 'ANTNTATTGTTAATGANTCGTTGTTTTACGTTT', 'REF': 'AGTGGATGGCGCCTGAGGCCTTGTTTGACCGAG', 'ID': '.'},
                   {'INFO': 'END=89690952;HOMLEN=3;HOMSEQ=TTA;SVLEN=5;SVTYPE=INS', 'FORMAT': 'GT:AD', 'CHROM': '10', 'POS': '89690952', 
                    'FILTER': 'PASS', 'QUAL': '.', 'READS': '1/1:11,67', 'ALT': 'TTTATC', 'REF': 'T', 'ID': '.'},
                   {'INFO': 'END=7579659;HOMLEN=25;HOMSEQ=CCCCAGCCCTCCAGGTCCCCAGCCC;SVLEN=-16;SVTYPE=DEL', 'FORMAT': 'GT:AD', 'CHROM': '17', 'POS': '7579643', 
                    'FILTER': 'PASS', 'QUAL': '.', 'READS': '0/1:35,26', 'ALT': 'C', 'REF': 'CCCCCAGCCCTCCAGGT', 'ID': '.'})
        self.refdata=(set([Interval(66763873, 66766605, {'bin': '136', 'exonFrames': '0,2,1,1,1,2,1,0,', 'exonNum': '1', 'cdsEndStat': 'cmpl', 'name2': 'AR', 'chrom': 'chrX', 'exonEnds': '66766604,66863249,66905968,66931531,66937464,66941805,66942826,66950461,', 'txEnd': '66950461', 'name': 'NM_000044', 'txStart': '66763873', 'cdsEnd': '66943683', 'score': '0', 'strand': '+', 'cdsStart': '66764988', 'cdsStartStat': 'cmpl', 'exonCount': '8', 'exonStarts': '66763873,66863097,66905851,66931243,66937319,66941674,66942668,66943527,'}), Interval(66766605, 66863096, {'bin': '136', 'exonFrames': '0,2,1,1,1,2,1,0,', 'cdsEndStat': 'cmpl', 'name2': 'AR', 'chrom': 'chrX', 'exonEnds': '66766604,66863249,66905968,66931531,66937464,66941805,66942826,66950461,', 'intronNum': '1', 'txEnd': '66950461', 'name': 'NM_000044', 'txStart': '66763873', 'cdsEnd': '66943683', 'score': '0', 'strand': '+', 'cdsStart': '66764988', 'cdsStartStat': 'cmpl', 'exonCount': '8', 'exonStarts': '66763873,66863097,66905851,66931243,66937319,66941674,66942668,66943527,'})]),
                     set([Interval(23634452, 23635328, {'bin': '765', 'exonFrames': '2,0,2,2,2,0,0,0,1,1,0,0,0,', 'cdsEndStat': 'cmpl', 'name2': 'PALB2', 'chrom': 'chr16', 'exonEnds': '23614990,23619333,23625412,23632799,23634451,23635415,23637718,23640596,23641790,23647655,23649273,23649450,23652678,', 'intronNum': '8', 'txEnd': '23652678', 'name': 'NM_024675', 'txStart': '23614482', 'cdsEnd': '23652478', 'score': '0', 'strand': '-', 'cdsStart': '23614779', 'cdsStartStat': 'cmpl', 'exonCount': '13', 'exonStarts': '23614482,23619184,23625324,23632682,23634289,23635329,23637556,23640524,23640960,23646182,23649170,23649390,23652430,'})]),
                      set([Interval(1807983, 1808054, {'bin': '598', 'exonFrames': '-1,0,1,1,1,0,1,0,1,0,2,1,1,0,0,2,2,0,', 'exonNum': '15', 'cdsEndStat': 'cmpl', 'name2': 'FGFR3', 'chrom': 'chr4', 'exonEnds': '1795192,1795770,1801250,1801539,1803263,1803470,1803752,1804791,1806247,1806696,1807203,1807396,1807667,1807900,1808054,1808410,1808661,1810599,', 'txEnd': '1810599', 'name': 'NM_001163213', 'txStart': '1795038', 'cdsEnd': '1808989', 'score': '0', 'strand': '+', 'cdsStart': '1795661', 'cdsStartStat': 'cmpl', 'exonCount': '18', 'exonStarts': '1795038,1795559,1800980,1801473,1803093,1803346,1803561,1804640,1806056,1806550,1807081,1807285,1807476,1807777,1807983,1808272,1808555,1808842,'}), Interval(1807983, 1808054, {'bin': '598', 'exonFrames': '-1,0,1,1,1,0,1,0,1,0,2,1,1,0,0,2,2,0,', 'exonNum': '15', 'cdsEndStat': 'cmpl', 'name2': 'FGFR3', 'chrom': 'chr4', 'exonEnds': '1795192,1795770,1801250,1801539,1803263,1803470,1803752,1805563,1806247,1806696,1807203,1807396,1807667,1807900,1808054,1808410,1808661,1810599,', 'txEnd': '1810599', 'name': 'NM_000142', 'txStart': '1795038', 'cdsEnd': '1808989', 'score': '0', 'strand': '+', 'cdsStart': '1795661', 'cdsStartStat': 'cmpl', 'exonCount': '18', 'exonStarts': '1795038,1795559,1800980,1801473,1803093,1803346,1803561,1805418,1806056,1806550,1807081,1807285,1807476,1807777,1807983,1808272,1808555,1808842,'}), Interval(1807983, 1808054, {'bin': '598', 'exonFrames': '-1,0,1,1,1,0,1,0,2,1,1,0,0,2,2,0,', 'exonNum': '13', 'cdsEndStat': 'cmpl', 'name2': 'FGFR3', 'chrom': 'chr4', 'exonEnds': '1795192,1795770,1801250,1801539,1803263,1803470,1803752,1806696,1807203,1807396,1807667,1807900,1808054,1808410,1808661,1810599,', 'txEnd': '1810599', 'name': 'NM_022965', 'txStart': '1795038', 'cdsEnd': '1808989', 'score': '0', 'strand': '+', 'cdsStart': '1795661', 'cdsStartStat': 'cmpl', 'exonCount': '16', 'exonStarts': '1795038,1795559,1800980,1801473,1803093,1803346,1803561,1806550,1807081,1807285,1807476,1807777,1807983,1808272,1808555,1808842,'})]))

    def testParseEvent(self):
        ''' Return the length and type of event, corrects end position if necessary '''
        #Test when start==stop but size > 1 (which only size >1 should ever hit this parser
        #Parse event returns negative value for DEL, but pindel_summary script converts this to abs for output.
        #Abs conversion is tested in the testPindelSummary test
        expected_output0=(-32, 'DEL',1808032)
        expected_output1=(5, 'INS', 89690953)
        expected_output2=(-16,'DEL',7579659)
        #test that RPL svtype becomes DEL
        self.assertEqual(pindel_summary.parse_event(self.data[0]), expected_output0)
        #Test that if start==stop, but size > 1, the stop is recalculated
        self.assertEqual(pindel_summary.parse_event(self.data[1]), expected_output1)
        #test that a 'normal' case processes correctly
        self.assertEqual(pindel_summary.parse_event(self.data[2]), expected_output2)
        

    def testDefineTranscripts(self):
        """Given the interval, set the gene, region and transcripts"""
    
        expected_output0=(['AR','AR'],['Exonic','Intronic'],['AR:NM_000044(exon 1)', 'AR:NM_000044(intron 1)'])
        expected_output1=(['PALB2'],['Intronic'],['PALB2:NM_024675(intron 8)'])
        expected_output2=(['FGFR3', 'FGFR3', 'FGFR3'],['Exonic','Exonic','Exonic'], ['FGFR3:NM_000142(exon 15)', 'FGFR3:NM_022965(exon 13)', 'FGFR3:NM_001163213(exon 15)'])
        #Test exonic
        self.assertEqual(pindel_summary.define_transcripts(self.refdata[0]), expected_output0)
        #Test intronic
        self.assertEqual(pindel_summary.define_transcripts(self.refdata[1]), expected_output1)
        # Test when start/stop cover multiple exons
        self.assertEqual(pindel_summary.define_transcripts(self.refdata[2]), expected_output2)


    def testPindelSummary(self):
        # Test when start/stop are in coding (ie normal case)
        # Test when start/stop are not incoding (ie intergenic case)
        # Test when start is in coding but stop isn't (edge case) 
        # Test when start is not in coding but stop is (edge case)
        # Test when a call covers exons and introns (edge case)
        simpletsv=os.path.join(pindel_testfiles, 'testing_output.txt')
        expected_output=os.path.join(pindel_testfiles, 'expected_output.txt')
        pindel_vcfs=[]
        for root, dirs, files in os.walk(pindel_testfiles):
            for file in files:
                if file.endswith(".vcf"):
                    pindel_vcfs.append(os.path.join(root, file))

        cmd=["./munge", "pindel_summary",self.refgene]+ [x for x in pindel_vcfs] +[ '-o',simpletsv ]
        subprocess.call(cmd)
        self.assertTrue(filecmp.cmp(expected_output, simpletsv))


