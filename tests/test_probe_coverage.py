"""
Test the subcommand scripts
"""
from os import path
import unittest
import logging
import subprocess
import filecmp

from __init__ import TestBase
import __init__ as config
log = logging.getLogger(__name__)


class Test_ProbeCoverage(TestBase):
    probe_testfiles = path.join(config.datadir, 'probe_coverage')    

    def test_probe_coverage(self):
        outdir = self.mkoutdir()
        tinput = path.join(self.probe_testfiles, 'test_input.tsv')
        toutput = path.join(self.probe_testfiles, "test_output.tsv")
        expected =  path.join(self.probe_testfiles, "expected_output.tsv")
        cmd=["./munge", "probe_coverage",  tinput, "-o" , toutput]
        subprocess.call(cmd)
        self.assertTrue(filecmp.cmp(toutput, expected))

