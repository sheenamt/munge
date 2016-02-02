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
from munging.subcommands import quality_metrics

from __init__ import TestBase
import __init__ as config
log = logging.getLogger(__name__)


class Test_QualityMetrics(TestBase):
    quality_testfiles = path.join(config.datadir, 'quality_metrics')    
    hs_file = path.join(quality_testfiles, '6037_E05_OPXv4_NA12878_HA0201.hs_metrics')    
    qm_file = path.join(quality_testfiles, '6037_E05_OPXv4_NA12878_HA0201.quality_metrics')       

    def test_quality(self):
        outdir = self.mkoutdir()
        simplecsv = path.join(outdir, "simple-qm.csv")
        expected =  path.join(self.quality_testfiles, "qm.csv")
        cmd=["munge", "quality_metrics", "-qm", self.qm_file, "-o" , simplecsv]
        subprocess.call(cmd)
        self.assertTrue(filecmp.cmp(simplecsv, expected))

    def test_hs(self):
        outdir = self.mkoutdir()
        simplecsv = path.join(outdir, "simple-hs.csv")
        expected =  path.join(self.quality_testfiles, "hs.csv")
        cmd=["munge", "quality_metrics", "-hs", self.hs_file, "-o" ,simplecsv]
        subprocess.call(cmd)
        self.assertTrue(filecmp.cmp(simplecsv, expected))

    def test_quality_and_hs(self):
        outdir = self.mkoutdir()
        simplecsv = path.join(outdir, "simple-hs-qm.csv")
        expected =  path.join(self.quality_testfiles, "hs-qm.csv")
        cmd=["munge", "quality_metrics", "-qm", self.qm_file,"-hs", self.hs_file, "-o" ,simplecsv]
        subprocess.call(cmd)
        self.assertTrue(filecmp.cmp(simplecsv, expected))

        
