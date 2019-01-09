"""
Test the breakdancer_summary script
"""

import subprocess
import filecmp
import logging
import os
from intervaltree import IntervalTree, Interval
from munging.subcommands import breakdancer_summary

from __init__ import TestBase
import __init__ as config
log = logging.getLogger(__name__)


breakdancer_testfiles = os.path.join(config.datadir, 'breakdancer')

class TestBreakdancerSummary(TestBase):
    """
    Test the breakdancer_summary script, which combines vcfs and annotates
    """
    def setUp(self):
        self.outdir = self.mkoutdir()
        self.refgene = os.path.join(breakdancer_testfiles, 'FGFR3_G6PD_IKBKG.txt')
        self.refdata={'chr4': IntervalTree([Interval(1795038, 1810599, {'bin': '598', 'exonEnds': '1795192,1795770,1801250,1801539,1803263,1803470,1803752,1806696,1807203,1807396,1807667,1807900,1808054,1808410,1808661,1810599,', 'exonFrames': '-1,0,1,1,1,0,1,0,2,1,1,0,0,2,2,0,', 'name': 'NM_022965', 'txStart': '1795038', 'exonCount': '16', 'cdsEndStat': 'cmpl', 'cdsEnd': '1808989', 'score': '0', 'name2': 'FGFR3', 'strand': '+', 'cdsStart': '1795661', 'cdsStartStat': 'cmpl', 'chrom': 'chr4', 'txEnd': '1810599', 'exonStarts': '1795038,1795559,1800980,1801473,1803093,1803346,1803561,1806550,1807081,1807285,1807476,1807777,1807983,1808272,1808555,1808842,'}),])}

    def testSetGeneEvent(self):
        """Given start position and gene interval tree, 
        return gene and event"""
        
        #Test normal
        expected_output0=('chr4:1806370','FGFR3')
        self.assertEqual(expected_output0, breakdancer_summary.set_gene_event('1806370', 'chr4', self.refdata))
        #Test intergenic
        expected_output1=('chr4:1795030','Intergenic')
        self.assertEqual(expected_output1, breakdancer_summary.set_gene_event('1795030', 'chr4', self.refdata))

    def testBreakdancerSummary(self):
        #Test that intergenic is processed correctly
        #Test that size filter is happening
        #Test that parsing is correct in general
        ctx=os.path.join(breakdancer_testfiles, 'BD.ctx')
        expected_output=os.path.join(breakdancer_testfiles, 'expected_output.txt')
        testing_output=os.path.join(breakdancer_testfiles, 'testing_output.txt')
        cmd=["munge", "breakdancer_summary", self.refgene, ctx, '-o', testing_output]
        subprocess.call(cmd)
        self.assertTrue(filecmp.cmp(expected_output, testing_output))
