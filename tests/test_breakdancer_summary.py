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
    Test the breakdancer_summary script, which parses the .ctx file, filters, and annotates
    """
    def setUp(self):
        self.outdir = self.mkoutdir()
        self.refgene = os.path.join(breakdancer_testfiles, 'BD_mini_refgene.txt')
        self.genes = os.path.join(breakdancer_testfiles, 'Breakdancer_genes.txt')
        self.ctx = os.path.join(breakdancer_testfiles, 'BD.ctx')
        
    def testUnfilteredBreakdancerSummary(self):
        #Test that intergenic is processed correctly
        #Test that size filter is happening
        #Test that parsing is correct in general
        expected_output=os.path.join(breakdancer_testfiles, 'expected_unfiltered_output.tsv')
        testing_output=os.path.join(self.outdir, 'testing_unfiltered_output.tsv')
        cmd=["./munge", "breakdancer_summary", self.refgene, self.ctx, '-o', testing_output]
        subprocess.call(cmd)
        self.assertTrue(filecmp.cmp(expected_output, testing_output))

    def testFilteredBreakdancerSummary(self):
        #Test that intergenic is processed correctly
        #Test that size filter is happening
        #Test that parsing is correct in general
        #Test that gene filter is happening properly
        expected_output=os.path.join(breakdancer_testfiles, 'expected_filtered_output.tsv')
        testing_output=os.path.join(self.outdir, 'testing_filtered_output.tsv')
        cmd=["./munge", "breakdancer_summary", self.refgene, self.ctx, '-g', self.genes, '-o', testing_output]
        subprocess.call(cmd)
        self.assertTrue(filecmp.cmp(expected_output, testing_output))