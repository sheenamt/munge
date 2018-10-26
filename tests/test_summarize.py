"""
Test the summarize_assay script
"""

import subprocess
import filecmp
import logging
import os

from munging.subcommands import summarize_assay

from __init__ import TestBase
import __init__ as config
log = logging.getLogger(__name__)


assay_testfiles = os.path.join(config.datadir, 'assay_files')


class TestSummaryAssay(TestBase):
    """
    Test the summarize_assay script, which gives overview of assay probes based on 
    UCSC refGene table info. 
    """
    def setUp(self):
        self.outdir = self.mkoutdir()

    def testSummarizeAssay(self):
        assay = os.path.join(assay_testfiles, 'test-bed-file.bed')
        pref_trans = os.path.join(assay_testfiles, 'test-gene-list.txt')
        refgene = os.path.join(assay_testfiles, 'test-refGene.bed')
        outdir = self.mkoutdir()
        bedtools='/mnt/disk2/com/container-images/bedtools-2.26.img'
        cmd=["munge", "summarize_assay", "--assay", assay, "--pref_trans", pref_trans, "--refgene", refgene, "--outdir", outdir, "--bedtools", bedtools]
        subprocess.call(cmd)

        #files made in this include "overall_summary.txt", "per_refseq_summary.txt", "merged_probes.bed"
        expected_overall=os.path.join(assay_testfiles, "expected-overall_summary.txt")
        expected_per_refseq=os.path.join(assay_testfiles, "expected-per_refseq_summary.txt")
        expected_merged_probes=os.path.join(assay_testfiles, "expected-merged_probes.bed")
        self.assertTrue(filecmp.cmp(expected_overall, os.path.join(outdir, "overall_summary.txt")))
        self.assertTrue(filecmp.cmp(expected_per_refseq, os.path.join(outdir, "per_refseq_summary.txt")))
        self.assertTrue(filecmp.cmp(expected_merged_probes, os.path.join(outdir, "merged_probes.bed")))
