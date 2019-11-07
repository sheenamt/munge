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
multiread_testfiles = os.path.join(config.datadir, 'pindel_multiread')

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
    def testParseEvent(self):
        ''' Return the length and type of event, corrects end position if necessary '''
        #Test when start==stop but size > 1 (which only size >1 should ever hit this parser
        #Parse event returns negative value for DEL, but pindel_summary script converts this to abs for output.
        #Abs conversion is tested in the testPindelSummary test
        expected_output0=(-32, 'DEL',1808033)
        expected_output1=(5, 'INS', 89690953)
        expected_output2=(-16,'DEL',7579660)
        #test that RPL svtype becomes DEL
        self.assertEqual(pindel_summary.parse_event(self.data[0]), expected_output0)
        #Test that if start==stop, but size > 1, the stop is recalculated
        self.assertEqual(pindel_summary.parse_event(self.data[1]), expected_output1)
        #test that a 'normal' case processes correctly
        self.assertEqual(pindel_summary.parse_event(self.data[2]), expected_output2)
        

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

    def testPindelSummaryMultiRead(self):
        # Test when start/stop are in coding (ie normal case)
        # Test when start/stop are not incoding (ie intergenic case)
        # Test when start is in coding but stop isn't (edge case) 
        # Test when start is not in coding but stop is (edge case)
        # Test when a call covers exons and introns (edge case)
        simpletsv=os.path.join(multiread_testfiles, 'multiread_testing_output.txt')
        expected_output=os.path.join(multiread_testfiles, 'multiread_expected_output.txt')
        pindel_vcfs=[]
        for root, dirs, files in os.walk(multiread_testfiles):
            for file in files:
                if file.endswith(".vcf"):
                    pindel_vcfs.append(os.path.join(root, file))

        cmd=["./munge", "pindel_summary", "--multi_reads", self.refgene]+ [x for x in pindel_vcfs] +[ '-o',simpletsv ]
        subprocess.call(cmd)
        self.assertTrue(filecmp.cmp(expected_output, simpletsv))


