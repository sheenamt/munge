"""
Test the breakdancer_summary script
"""

import subprocess
import filecmp
import logging
import os

from munging.subcommands import breakdancer_summary

from __init__ import TestBase
import __init__ as config
log = logging.getLogger(__name__)


refgene_testfiles = os.path.join(config.datadir, 'refseq_filter')
breakdancer_testfiles = os.path.join(config.datadir, 'breakdancer')

class TestBreakdancerSummary(TestBase):
    """
    Test the breakdancer_summary script, which combines vcfs and annotates
    """
    def setUp(self):
        self.outdir = self.mkoutdir()
        self.refgene = os.path.join(refgene_testfiles, 'expected_refgene.txt')

    def testBreakdancerSummary(self):
        #Test that if a gene has a preferred transcript, that one is chosen
        #Test that if a gene doesn't have a preferred transcript, one line of data is still output
        #Test that overlapping genes are allowed 
        breakdancer_vcfs=[]
        for root, dirs, files in os.walk(breakdancer_testfiles):
            for file in files:
                if file.endswith(".vcf"):
                    breakdancer_vcfs.append(os.path.join(root, file))
        breakdancer_vcfs=str(breakdancer_testfiles)+'/*vcf'
        cmd=["munge", "breakdancer_summary", self.refgene,breakdancer_vcfs  ] #, '-o',simplecsv ]
        subprocess.call(cmd)
#        self.assertTrue(filecmp.cmp(self.expected_output, simplecsv))
