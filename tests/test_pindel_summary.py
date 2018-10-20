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

    def testParseEvent(self):
    ''' Return the length and type of event, corrects end position if necessary '''
    #Test when start==stop but size > 1 (which only size >1 should ever hit this parser
    pass

    def testDefineTranscripts(self):
    """Given the interval, set the gene, region and transcripts"""
    pass
    # Test when start/stop cover multiple exons


    def testPindelSummary(self):
        # Test when start/stop are in coding
        # Test when start/stop are not incoding
        # Test when start is in coding but stop isn't
        # Test when start is not in coding but stop is

        pindel_vcfs=[]
        for root, dirs, files in os.walk(pindel_testfiles):
            for file in files:
                if file.endswith(".vcf"):
                    pindel_vcfs.append(os.path.join(root, file))
        pindel_vcfs=str(pindel_testfiles)+'/*vcf'
        cmd=["munge", "pindel_summary", self.refgene,pindel_vcfs  ] #, '-o',simplecsv ]
        subprocess.call(cmd)
#        self.assertTrue(filecmp.cmp(self.expected_output, simplecsv))


