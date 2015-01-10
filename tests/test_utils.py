"""
Test the utils functions
"""

import os
from os import path
import unittest
import logging
import pprint
import sys
import json

from munging.utils import munge_path, munge_pfx

from __init__ import TestBase
import __init__ as config
log = logging.getLogger(__name__)

class TestUtils(TestBase):

    def setUp(self):
        self.run=munge_path('/home/genetics/data/2013-03-21_HiSeq_OncoPlex_Run14_v2')
        self.outdir = self.mkoutdir()

    def tearDown(self):
        pass

    def testMungePathSplitting(self):
        """
        Tests correct splitting of path into assay, run, machine, capture
        """
        self.assertEquals(self.run['assay'], 'oncoplex')
        self.assertEquals(self.run['run'], '2013-03-21_run14')
        self.assertEquals(self.run['machine'], 'hiseq')
        self.assertEquals(self.run['capture'], 'v2')


    def testMungePathLowerCase(self):
        """
        Tests correct lowercase of path into assay, run, machine
        """
        self.assertNotEqual(self.run['assay'], 'OncoPlex')
        self.assertNotEqual(self.run['run'], '2013-03-21_Run14')
        self.assertNotEqual(self.run['machine'], 'HiSeq')

    def testMungeWrongPath(self):
        """
        Test that a ValueError is raised if wrong path is given
        """
        run=munge_path('/home/genetics/data/Broca_v6')
        self.assertRaises(ValueError)

    def testMungePFX(self):
        real_info={'control': 'NA12878', 
                   'machine-run': 'HA0201', 
                   'library-version': 'OPXv4', 
                   'well': 'E05', 
                   'run': '60', 
                   'sample_id': '6037', 
                   'pfx': '6037_E05_NA12878_OPXv4', 
                   'mini-pfx': '6037_NA12878'}
        
        test_info=munge_pfx('6037_E05_OPXv4_NA12878_HA0201')
        self.assertDictEqual(real_info, test_info)
    def testCheckControl(self):
        pass    