"""
Test the subcommand scripts
"""
import os
from os import path
import unittest
import logging
import pprint
import csv
import sys
import json
import subprocess
import filecmp
from munging.subcommands import parse_monoseq

from __init__ import TestBase
import __init__ as config
log = logging.getLogger(__name__)
                
class TestMonoSeqParser(TestBase):
    """
    Test the monoseq parser, which makes the ugly format human readable
    """
    analysis_testfiles = os.path.join(config.datadir, '101010_HA0000_OncoPlex1','output','CFDNA')
    manifest = os.path.join(config.datadir, '101010_HA0000_OncoPlex1', 'configs','pipeline-manifest-CFDNA.csv')
    monoseq_testfiles = os.path.join(config.datadir, '101010_HA0000_OncoPlex1','data')

    def testMonoSeq(self):
        sample_id = 'A-NA12878-LMG-240-GEN0110-GLTv1'
        monoseq_input = path.join(self.analysis_testfiles, '{}.CFDNA.monoseq').format(sample_id)
        monoseq_output = path.join(self.analysis_testfiles, '{}.CFDNA.PolyHunter_Analysis.txt').format(sample_id)
        outdir = self.mkoutdir()
        simplecsv = path.join(outdir, "simple-monoseq.csv")
        expected =  monoseq_output
        cmd=["./munge", "parse_monoseq", monoseq_input, path.join(self.monoseq_testfiles, 'CFTR.xml'), simplecsv]
        subprocess.call(cmd)
        self.assertTrue(filecmp.cmp(simplecsv, expected))

    def testMonoSeqTopLevel(self):
        outdir = self.mkoutdir()
        simplecsv = path.join(outdir, "simple-monoseq-top-level.csv")
        monoseq_output = path.join(self.analysis_testfiles, 'GEN0110-GLTv1.Combined_PolyHunter.txt')
        expected =  monoseq_output
        cmd=["./munge", "parse_monoseq_top_level", self.analysis_testfiles, self.manifest, '-o', simplecsv]
        subprocess.call(cmd)
        self.assertTrue(filecmp.cmp(simplecsv, expected))
