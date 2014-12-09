"""
Test the annotation functions
"""

import os
from os import path
import unittest
import logging
import pprint
import sys
import json

from munging.annotation import get_location
from munging.annotation import multi_split
from munging.annotation import split_string_in_two
from munging.annotation import split_chr_loc
from munging.annotation import build_variant_id
from munging.annotation import fix_pfx

from __init__ import TestBase
import __init__ as config
log = logging.getLogger(__name__)

class TestAnnotation(TestBase):

    def setUp(self):
        self.outdir = self.mkoutdir()

    def tearDown(self):
        pass

    def testGetLocation01(self):
        """
        Tests string arguments
        """
        loc = get_location('1', '1', '100')
        self.assertEquals(loc, 'chr1:1-100')
        loc = get_location('1', '1', '1')
        self.assertEquals(loc, 'chr1:1')

    def testGetLocation02(self):
        """
        Tests integer arguments
        """
        loc = get_location(1, 1, 100)
        self.assertEquals(loc, 'chr1:1-100')
        loc = get_location(1, 1, 1)
        self.assertEquals(loc, 'chr1:1')

    def testSplitString(self):
        """
        Tests spliting a string given string of split characters
        Used to split path names in crawlers
        """
        result=multi_split('/home/genetics/data/run_info', '/_')
        self.assertEquals(result, ['home','genetics','data','run','info'])


    def testSplitChrLoc(self):
        """
        Tests spliting the chr_loc from chr1:1-100 to 1, 1, 100
        Returns strings, not integers
        """
        result01=split_chr_loc('chr1:1-100')
        result02=split_chr_loc('1:1-100')
        self.assertEquals(result01, ('1','1','100'))
        self.assertEquals(result02, ('1','1','100'))

    def testSplitStringinTwo(self):
        """
        Tests spliting a string into two on | or ,
        Return -1, -1 if input is None
        Used to split columns of data in summary
        """
        result01=split_string_in_two("1|99")
        result02=split_string_in_two(None)
        self.assertEquals(result01[0], '1')
        self.assertEquals(result01[1], '99')
        self.assertEquals(result02[0], '-1')
        self.assertEquals(result02[1], '-1')

    def testBuildVariantID(self):
        """
        Tests creation of variant id from
        chrX:start-stop, ref_base, var_base or
        chrX:start, ref_base, var_base
        """
        data1=['chrX:12321', 'A', 'T']
        data2=['chrX:1234-1256', 'G', 'C']
        result01=build_variant_id(data1)
        result02=build_variant_id(data2)
        self.assertEquals(result01, 'X_12321_12321_A_T')
        self.assertEquals(result02, 'X_1234_1256_G_C')

    def testFixPfx(self):
        self.assertEqual(fix_pfx('48_A03_BROv7-HA0186-NA12878'), '48_A03_BROv7_HA0186_NA12878')
        self.assertEqual(fix_pfx('48_A03_BROv7-HA0186-NA12878 '), '48_A03_BROv7_HA0186_NA12878')
        self.assertEqual(fix_pfx('UNK-124-455'), 'UNK_124_455')
        self.assertEqual(fix_pfx('LMG-240'), 'LMG240')

