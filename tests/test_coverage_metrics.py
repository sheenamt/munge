'''
Test the coverage_metrics subcommand
'''
import subprocess
import filecmp
import StringIO
import sys
import os
import unittest
import logging
import pandas as pd
from munging.subcommands import coverage_metrics

from __init__ import TestBase
import __init__ as config
log = logging.getLogger(__name__)

cov_testfiles = os.path.join(config.datadir, 'coverage')
class TestCoverageMetrics(TestBase):
    '''
    Test the coverage metrics script
    '''
    
    def setUp(self):
        self.outdir = self.mkoutdir()
        self.target_cov_input=os.path.join(cov_testfiles,'BRCA1_PerBaseCoverage.txt')
        self.refgene_input=os.path.join(cov_testfiles,'test-refGene.txt')

    def testParsePerBase50(self):
        """Test parsing of per base output from bedtools with a failing coverage threshold of 50"""
        expected=os.path.join(cov_testfiles,'expected_output_t50.tsv')
        simpletsv=os.path.join(self.outdir, "simple-brca1-t50.tsv")
        cmd=["./munge", "coverage_metrics", self.target_cov_input, self.refgene_input, "-o", simpletsv]
        subprocess.call(cmd)
        self.assertTrue(filecmp.cmp(simpletsv, expected))

    def testParsePerBase100(self):
        """Test parsing of per base output from bedtools with a failing coverage threshold of 100"""
        expected=os.path.join(cov_testfiles,'expected_output_t100.tsv')
        simpletsv=os.path.join(self.outdir, "simple-brca1-t100.tsv")
        cmd=["./munge", "coverage_metrics", self.target_cov_input, self.refgene_input, "-t", "100", "-o", simpletsv]
        subprocess.call(cmd)
        self.assertTrue(filecmp.cmp(simpletsv, expected))


