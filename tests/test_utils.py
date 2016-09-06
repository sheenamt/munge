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

    def testMungePFX1(self):
        real_info={'control': 'NA12878', 
                   'machine-run': 'HA0201', 
                   'library-version': 'OPXv4', 
                   'well': 'E05', 
                   'run': '60', 
                   'sample_id': '6037', 
                   'pfx': '6037_E05_OPXv4_NA12878_HA0201',
                   'assay':'OncoPlex',
                   'mini-pfx': '6037_NA12878'}
        
        test_info=munge_pfx('6037_E05_OPXv4_NA12878_HA0201')
        self.assertDictEqual(real_info, test_info)
    def testMungePFX2(self):
        real_info={'sample_id':'LMG-240',
                   'pfx': 'LMG-240',
                   'assay':'Coloseq',
                   'mini-pfx': 'LMG-240'}
        test_info=munge_pfx('LMG-240')
        self.assertDictEqual(real_info, test_info)

    def testMungePFX3(self):
        real_info={'sample_id': '6037', 
                   'pfx': '6037',
                   'assay':'Coloseq',
                   'well':'01',
                   'mini-pfx': '6037'}
        test_info=munge_pfx('6037_01_BROv8')

    def testMungePFX4(self):
        real_info={'sample_id': '6037', 
                   'pfx': '6037',
                   'assay':'MSI-PLUS',
                   'mini-pfx': '6037'}
        test_info=munge_pfx('6037_MSI-Plus')

    def testMungePFX5(self):
        real_info={'sample_id': '6037', 
                   'pfx': '6037',
                   'assay':'hotspot-hereditary',
                   'mini-pfx': '6037'}
        test_info=munge_pfx('6037_GLT06')

    def testMungePFX6(self):
        real_info={'sample_id': '6037', 
                   'pfx': '6037',
                   'assay':'hotspot-heme',
                   'mini-pfx': '6037'}
        test_info=munge_pfx('6037_HP06')

    def testMungePath(self):
        test_id1=munge_path('140915_HA000_ColoTestFiles')
        test_id2=munge_path('testfiles/140915_MA0001_OncoPlexKapa')
        test_id3=munge_path('testfiles/160727_MB0745_GEN005-GLTv1')
        test_id4=munge_path('testfiles/160820_MD0020_HP66-HHv1')
        test_id5=munge_path('testfiles/160209_MA0829_GEN06-STHv1')

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

        self.assertEqual(test_id3,({'date':'2016-07-27',
                                    'run': 'MB0745', 
                                    'machine': 'miseq',
                                    'assay':'hotspot-hereditary', 
                                    'prep_type':'truseq',
                                    'project': 'gen005-gltv1'}))

        self.assertEqual(test_id4,({'date':'2016-08-20',
                                    'run': 'MD0020', 
                                    'machine': 'miseq',
                                    'assay':'hotspot-heme', 
                                    'prep_type':'truseq',
                                    'project': 'hp66-hhv1'}))

        self.assertEqual(test_id5,({'date':'2016-02-09',
                                    'run': 'MA0829', 
                                    'machine': 'miseq',
                                    'assay':'hotspot-heme', 
                                    'prep_type':'truseq',
                                    'project': 'gen06-sthv1'}))

    def testMungedate(self):
        test_id1=munge_date('140915')
        test_id2=munge_date('2014-09-09')

        self.assertEqual(test_id1,'2014-09-15')
        self.assertEqual(test_id2,'2014-09-09')
