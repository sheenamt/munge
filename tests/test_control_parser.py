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

from munging.subcommands import control_parser

from __init__ import TestBase
import __init__ as config
log = logging.getLogger(__name__)

control_testfiles = path.join(config.datadir, 'control_parser')

class TestControlParser(TestBase):
    """
    Test the control_parser which is used to check the control sample against
    SNVs found in NA12878 from Complete Genomics, 1000G and our samples
    """

    def testControlParserMatch(self):
        """
        Test when the hapmap vcf output and the pipeline output are 100% in agreement
        """
        controlfname = path.join(control_testfiles, 'HAPMAP_vcf_variant_function')
        runfname = path.join(control_testfiles, 'HAPMAP_pipeline_output_full')
        outdir = self.mkoutdir()
        simplecsv = path.join(outdir, "hapmap_qc_full.csv")
        expected =  path.join(control_testfiles, "hapmap_expected_full.csv")
        cmd=["munge", "control_parser", controlfname, runfname, "-o", simplecsv]
        subprocess.call(cmd)
        self.assertTrue(filecmp.cmp(simplecsv, expected))


    def testControlParserMatch2(self):
        """
        Test when the hapmap vcf output has more output than the pipeline data
        """
        controlfname = path.join(control_testfiles, 'HAPMAP_vcf_variant_function')
        runfname = path.join(control_testfiles, 'HAPMAP_pipeline_output_missing')
        outdir = self.mkoutdir()
        simplecsv = path.join(outdir, "hapmap_qc_missing.csv")
        expected =  path.join(control_testfiles, "hapmap_expected_missing.csv")
        cmd=["munge", "control_parser", controlfname, runfname, "-o", simplecsv]
        subprocess.call(cmd)
        self.assertTrue(filecmp.cmp(simplecsv, expected))


