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

    def testParsePerBase1(self):
        '''Test parsing of per base output from bedtools
        '''
        expected=os.path.join(cov_testfiles,'expected_output.txt')
        simplecsv=os.path.join(self.outdir, "simple-brca1-100.csv")
        cmd=["./munge", "coverage_metrics", self.target_cov_input, self.refgene_input, "-o", simplecsv]
        subprocess.call(cmd)
        self.assertTrue(filecmp.cmp(simplecsv, expected))


