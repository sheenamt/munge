"""
Test the pindel_summary script
"""

import subprocess
import filecmp
import logging
import os

from munging.subcommands import pindel_summary

from __init__ import TestBase
import __init__ as config
log = logging.getLogger(__name__)


refgene_testfiles = os.path.join(config.datadir, 'refseq_filter')
pindel_testfiles = os.path.join(config.datadir, 'pindel')

class TestPindelSummary(TestBase):
    """
    Test the pindel_summary script, which combines vcfs and annotates
    """
    def setUp(self):
        self.outdir = self.mkoutdir()
        self.refgene = os.path.join(refgene_testfiles, 'expected_refgene.txt')

    def testPindelSummary(self):
        #Test that if a gene has a preferred transcript, that one is chosen
        #Test that if a gene doesn't have a preferred transcript, one line of data is still output
        #Test that overlapping genes are allowed 
        pindel_vcfs=[]
        for root, dirs, files in os.walk(pindel_testfiles):
            for file in files:
                if file.endswith(".vcf"):
                    pindel_vcfs.append(os.path.join(root, file))
        pindel_vcfs=str(pindel_testfiles)+'/*vcf'
        cmd=["munge", "pindel_summary", self.refgene,pindel_vcfs  ] #, '-o',simplecsv ]
        subprocess.call(cmd)
#        self.assertTrue(filecmp.cmp(self.expected_output, simplecsv))
