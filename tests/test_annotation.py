"""
Test the annotation functions
"""

import logging
import os
import __init__ as config
from intervaltree import Interval, IntervalTree

from munging.annotation import get_location
from munging.annotation import multi_split
from munging.annotation import split_chr_loc
from munging.annotation import split_string_in_two
from munging.annotation import build_variant_id
from munging.annotation import pfx_ok
from munging.annotation import fix_pfx
from munging.annotation import _fix
from munging.annotation import GenomeIntervalTree
from munging.annotation import UCSCTable
from munging.annotation import IntervalMakers


from __init__ import TestBase
 
log = logging.getLogger(__name__)

class TestAnnotation(TestBase):

    def setUp(self):
        self.outdir = self.mkoutdir()
        self.refgene = os.path.join(config.datadir, 'pindel', 'refgene_test.txt')

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

    def test_fix(self):
        ok_interval=Interval(1,10,('1-10'))
        bkwds_interval=Interval(10,1,('10-1'))
        
        correct_interval_ok=Interval(1,10,('1-10'))
        correct_interval_bkwds=Interval(10,11,('10-1'))
        self.assertEqual(_fix(ok_interval), correct_interval_ok)
        self.assertEqual(_fix(bkwds_interval), correct_interval_bkwds)

    def testGenomeIntervalTree(self):
        data=IntervalTree()
        exons=GenomeIntervalTree.from_table(open(self.refgene, 'r'), parser=UCSCTable.REF_GENE, mode='exons')

        #642     NM_000546       chr17   -       7571719 7590868 7572926 7579912 11      7571719,7573926,7576852,7577018,7577498,7578176,7578370,7579311,7579699,7579838,7590694,        7573008,7574033,7576926,7577155,7577608,7578289,7578554,7579590,7579721,7579940,7590868,
        for start,end,data in exons['chr17'].search(int(7577100)):
            if data['name']=='NM_000546':
                rev_exon=data
        for start,end,data in exons['chr17'].search(int(7577700)):
            if data['name']=='NM_000546':
                rev_intron=data

        #TP53:NM_000546, reverse strand, position 7577100 is in exon8,position 7577100 is in intron6
        self.assertEqual(rev_exon['exonNum'],'8')
        self.assertEqual(rev_intron['intronNum'],'6')

        ##3' UTR is  7590868-7579912
        for start,end,data in exons['chr17'].search(int(7571719)):
            if data['name']=='NM_000546':
                rev_3_utr_start=data

        for start,end,data in exons['chr17'].search(int(7572926)):
            if data['name']=='NM_000546':
                rev_3_utr_stop=data

        #5' UTR is 7571719-7572926 
        for start,end,data in exons['chr17'].search(int(7590868)):
            if data['name']=='NM_000546':
                rev_5_utr_start=data

        for start,end,data in exons['chr17'].search(int(7579912)):
            if data['name']=='NM_000546':
                rev_5_utr_stop=data
 
        #TP53:NM_000546, reverse strand, 3' UTR is 7571719-7572926, 5' UTR is 7590868-7579912
        self.assertEqual(rev_3_utr_start['UTR'],'3')
        self.assertEqual(rev_3_utr_stop['UTR'],'3')
        self.assertEqual(rev_5_utr_start['UTR'],'5')
        self.assertEqual(rev_5_utr_stop['UTR'],'5')
        


        for start,end,data in exons['chrX'].search(int(66763880)):
            if data['name']=='NM_000044':
                forward_exon=data
        for start,end,data in exons['chrX'].search(int(66942830)):
            if data['name']=='NM_000044':
                forward_intron=data
        #AR:NM_000044, forward strand, position 66763880  is in exon1,position 66942830 is in intron7
        self.assertEqual(forward_exon['exonNum'],'1')
        self.assertEqual(forward_intron['intronNum'],'7')

        
        #3' UTR is 66950461-66943683
        for start,end,data in exons['chrX'].search(int(66950461)):
            if data['name']=='NM_000044':
                utr_3_start=data

        #Test 3' UTR
        for start,end,data in exons['chrX'].search(int(66943683)):
            if data['name']=='NM_000044':
                utr_3_stop=data

        #5' UTR - 66763873-66764988
        for start,end,data in exons['chrX'].search(int(66763873)):
            if data['name']=='NM_000044':
                utr_5_start=data

        # Test 5' UTR
        for start,end,data in exons['chrX'].search(int(66764988)):
            if data['name']=='NM_000044':
                utr_5_stop=data

        self.assertEqual(utr_3_start['UTR'],'3')
        self.assertEqual(utr_3_stop['UTR'],'3')
        self.assertEqual(utr_5_start['UTR'],'5')
        self.assertEqual(utr_5_stop['UTR'],'5')
  
        # #Test 0 based vs 1 based entry for start and stop
        # # Test first exon and last exon to make sure we are calling correctly the difference between UTR and exon
        # #1-10 UTR, 11-20 is exon1:
        # #10 == UTR, 11=exon1, 20=exon1, 21=intron1
