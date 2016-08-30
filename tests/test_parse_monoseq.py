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
    analysis_testfiles = path.join(config.datadir, 'analysis_files')    
    monoseq_testfiles = path.join(config.datadir, 'monoseq')    
    monoseq_input = path.join(monoseq_testfiles, 'NA12785A-GLT005.CFDNA.monoseq')
    monoseq_output = path.join(analysis_testfiles, 'NA12785A-GLT005.CFDNA.PolyHunter_Analysis.txt')

    def testMonoSeq(self):
        outdir = self.mkoutdir()
        simplecsv = path.join(outdir, "simple-monoseq.csv")
        expected =  self.monoseq_output
        cmd=["munge", "parse_monoseq", self.monoseq_input, path.join(self.monoseq_testfiles, 'CFTR-polyT.xml'), simplecsv]
        subprocess.call(cmd)
        self.assertTrue(filecmp.cmp(simplecsv, expected))
