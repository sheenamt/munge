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

from munging.utils import munge_path

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
