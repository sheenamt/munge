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

from munging.utils import munge_path, munge_pfx, munge_date

from __init__ import TestBase
import __init__ as config
log = logging.getLogger(__name__)

class TestUtils(TestBase):

    def setUp(self):
        self.run=munge_path('/home/genetics/data/130321_HA00211_OncoPlex64')
        self.outdir = self.mkoutdir()

    def testMungePFX(self):
        real_info={'control': 'NA12878', 
                   'machine-run': 'HA0201', 
                   'library-version': 'OPXv4', 
                   'well': 'E05', 
                   'run': '60', 
                   'sample_id': '6037', 
                   'pfx': '6037_E05_OPXv4_NA12878_HA0201',
                   'assay':'oncoplex',
                   'mini-pfx': '6037_NA12878'}
        
        test_info=munge_pfx('6037_E05_OPXv4_NA12878_HA0201')
        self.assertDictEqual(real_info, test_info)

    def testCheckControl(self):
        pass    


    def testMungePath(self):
        test_id1=munge_path('140915_HA000_ColoTestFiles')
        test_id2=munge_path('testfiles/140915_MA0001_OncoPlexKapa')

        self.assertEqual(test_id1,({'date':'2014-09-15',
                                    'run': 'HA000', 
                                    'machine': 'hiseq',
                                    'assay':'coloseq',
                                    'project': 'colotestfiles',
                                    'prep_type':'sure_select'}))

        self.assertEqual(test_id2,({'date':'2014-09-15',
                                    'run': 'MA0001', 
                                    'machine': 'miseq',
                                    'assay':'oncoplex', 
                                    'prep_type':'kapa',
                                    'project': 'oncoplexkapa'}))

    def testMungedate(self):
        test_id1=munge_date('140915')
        test_id2=munge_date('2014-09-09')

        self.assertEqual(test_id1,'2014-09-15')
        self.assertEqual(test_id2,'2014-09-09')
