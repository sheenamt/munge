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

from munging.utils import munge_path, munge_pfx, munge_date, munge_samples

from __init__ import TestBase
import __init__ as config
log = logging.getLogger(__name__)

class TestUtils(TestBase):

    def setUp(self):
        self.run='/home/genetics/data/130321_HA00211_OncoPlex64'
        self.outdir = self.mkoutdir()

    def testMungePFX1(self):
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
    def testMungePFX2(self):
        real_info={'sample_id':'LMG-240',
                   'pfx': 'LMG-240',
                   'assay':'coloseq',
                   'mini-pfx': 'LMG-240'}
        test_info=munge_pfx('LMG-240')
        self.assertDictEqual(real_info, test_info)

    def testMungePFX3(self):
        real_info={'sample_id': '6037', 
                   'pfx': '6037',
                   'assay':'coloseq',
                   'well':'01',
                   'mini-pfx': '6037'}
        test_info=munge_pfx('6037_01_BROv8')

    def testMungePFX4(self):
        real_info={'sample_id': '6037', 
                   'pfx': '6037',
                   'assay':'msi-plus',
                   'mini-pfx': '6037'}
        test_info=munge_pfx('6037_MSI-Plus')

    def testMungePath(self):
        test_id1=munge_path('140915_HA000_ColoTestFiles')
        test_id2=munge_path('testfiles/140915_MA0001_OncoPlexKapa')
        test_id3=munge_path('testfiles/140915_MA0001_OncoPlexKapa/output')

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

        self.assertEqual(test_id3,({'date':'2014-09-15',
                                    'run': 'MA0001', 
                                    'machine': 'miseq',
                                    'assay':'oncoplex', 
                                    'prep_type':'kapa',
                                    'project': 'oncoplexkapa'}))

        self.assertRaises(ValueError,munge_path,'testfiles/140915_MA0001_MSIPlus')
        self.assertRaises(ValueError,munge_path,'testfiles/140915_NPM1_0123')
        self.assertRaises(ValueError,munge_path,'testfiles/1405_07')
        self.assertRaises(ValueError,munge_path,'150813_MA0089_MSIplus_150811')
        self.assertRaises(ValueError,munge_path,'08-19-15_NPM1NG2')
        self.assertRaises(ValueError,munge_path,'150728_NPM1_PP7')
        self.assertRaises(ValueError,munge_path,'TruSeq06')
        self.assertRaises(ValueError,munge_path,'2015_MSI-PLUS_Validation')

    def testMungedate(self):
        test_id1=munge_date('140915')
        test_id2=munge_date('2014-09-09')

        self.assertEqual(test_id1,'2014-09-15')
        self.assertEqual(test_id2,'2014-09-09')


class TestManifest(TestBase):
    def testMungeSamples(self):
        pth1='testfiles/'
        for sample in munge_samples(pth1):
            if sample['pfx']=='7-LMG240':
                self.assertTrue(sample['is_control'])
                self.assertEqual(sample['run'], '7')
                self.assertEqual(sample['project'], 'testfiles')
            if sample['pfx']=='7-B1':
                self.assertFalse(sample['is_control'])
                self.assertEqual(sample['run'], '7')
                self.assertEqual(sample['project'], 'testfiles')
                
                
