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
from munging.annotation import get_exons


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
        result01=split_chr_loc('chr1:1-100 ')
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
        Tests creation of variant id and read count from
        chrX:start-stop, ref_base, var_base or
        chrX:start, ref_base, var_base
        """
        #A bunch of empty 'columns' to represent the data in the SNP tab
        data1=['chrX:12321', 'A', 'T','','','','','','','','','','','','','10','11']
        data2=['chrX:1234-1256', 'G', 'C','','','','','','','','','','','','','12','13']
        result01,ref_reads01,varreads01=build_variant_id(data1)
        result02,ref_reads02,varreads02=build_variant_id(data2)
        self.assertListEqual([result01,ref_reads01,varreads01],['X_12321_12321_A_T','10','11'])
        self.assertListEqual([result02,ref_reads02,varreads02], ['X_1234_1256_G_C','12','13'])

        

    def testFixPfx(self):
        self.assertEqual(fix_pfx('48_A03_BROv7-HA0186-NA12878'), '48_A03_BROv7_HA0186_NA12878')
        self.assertEqual(fix_pfx('48_A03_BROv7-HA0186-NA12878 '), '48_A03_BROv7_HA0186_NA12878')
        self.assertEqual(fix_pfx('UNK-124-455'), 'UNK_124_455')
        self.assertEqual(fix_pfx('LMG-240'), 'LMG240')

    def testGetExons(self):
        """Test the creation of exon/intron partion, 
        which maps exon/intron number to coordinates"""

        starts='66763873,66863097,66905851,66931243,66937319,66941674,66942668,66943527,'
        ends='66766604,66863249,66905968,66931531,66937464,66941805,66942826,66950461,'

        forward_expected=[('ex1', 66763873, 66766604), ('ex2', 66863097, 66863249), ('ex3', 66905851, 66905968), ('ex4', 66931243, 66931531), ('ex5', 66937319, 66937464), ('ex6', 66941674, 66941805), ('ex7', 66942668, 66942826), ('ex8', 66943527, 66950461), ('int1', 66766605, 66863096), ('int2', 66863250, 66905850), ('int3', 66905969, 66931242), ('int4', 66931532, 66937318), ('int5', 66937465, 66941673), ('int6', 66941806, 66942667), ('int7', 66942827, 66943526)]
        reverse_expected= [('ex1', 66943527, 66950461), ('ex2', 66942668, 66942826), ('ex3', 66941674, 66941805), ('ex4', 66937319, 66937464), ('ex5', 66931243, 66931531), ('ex6', 66905851, 66905968), ('ex7', 66863097, 66863249), ('ex8', 66763873, 66766604), ('int1', 66942827, 66943526), ('int2', 66941806, 66942667), ('int3', 66937465, 66941673), ('int4', 66931532, 66937318), ('int5', 66905969, 66931242), ('int6', 66863250, 66905850), ('int7', 66766605, 66863096)]
        forward_output=get_exons(starts, ends, '+')
        reverse_output=get_exons(starts, ends, '-')
        self.assertEqual(forward_expected, forward_output)
        self.assertEqual(reverse_expected, reverse_output)

