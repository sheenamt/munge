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
from munging.annotation import define_transcripts


from __init__ import TestBase
 
log = logging.getLogger(__name__)

class TestAnnotation(TestBase):

    def setUp(self):
        self.outdir = self.mkoutdir()
        self.refgene = os.path.join(config.datadir, 'pindel', 'refgene_test.txt')
        self.exons=GenomeIntervalTree.from_table(open(self.refgene, 'r'), parser=UCSCTable.REF_GENE, mode='exons')

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

    def testDefineTranscripts(self):
        """Given the interval, set the gene, region and transcripts"""
        #Test Exonic region, when only exonic (66905851 is exon3 start, 1based, 66905968 is exon3 end 1based )
        ar_exon3=self.exons['chrX'].search(int(66905851),int(66905968))
        expected0=(['AR'], ['EXONIC'], ['AR:NM_000044(exon 03)', 'AR:NM_001011645(exon 03)'])
        self.assertEqual(define_transcripts(ar_exon3),expected0)


        #Test Exonic region, when exonic and intronic (exon2 + intron3)
        ar_exon2_intronic3=self.exons['chrX'].search(int(66863100), int(66905970))
        expected1=(['AR'], ['EXONIC'], ['AR:NM_000044(exon 02 - intron 03)', 'AR:NM_001011645(exon 02 - intron 03)'])
        self.assertEqual(define_transcripts(ar_exon2_intronic3), expected1)
        #exon 2, intron 3

        #Test Exonic region, when exonic and intronic and UTR
        ar_utr_exon2_intronic3=self.exons['chrX'].search(int(66764987), int(66905970))
        expected3=(['AR'], ['EXONIC'], ['AR:NM_000044(UTR - intron 03)', 'AR:NM_001011645(UTR - intron 03)'])
        self.assertEqual(define_transcripts(ar_utr_exon2_intronic3), expected3)

        ##Test intronic
        ar_intronic=self.exons['chrX'].search(int(66905970), int(66905975))
        expected4=(['AR'], ['INTRONIC'], ['AR:NM_000044(intron 03)', 'AR:NM_001011645(intron 03)'])
        self.assertEqual(define_transcripts(ar_intronic), expected4)

        #Test 5 UTR on - 
        #txStart (532241) and cdsStart (532635)
        hras_utr=self.exons['chr11'].search(int(532241), int(532635))
        expected5=(['HRAS'],['UTR'],['HRAS:NM_005343(UTR)'])
        self.assertEqual(define_transcripts(hras_utr),expected5)
        
        #Test 3 UTR on - 
        #cdsEnd(534322) and txEnd (535567) 
        hras_3_utr=self.exons['chr11'].search(int(534323), int(535568))
        expected6=(['HRAS'],['UTR'],['HRAS:NM_005343(UTR)'])
        self.assertEqual(define_transcripts(hras_3_utr),expected6)

        #Test UTR on +
        #txStart (66763873) and cdsStart (66764988)
        ar_utr=self.exons['chrX'].search(int(66763873), int(66764988))
        expected7=(['AR'],['UTR'],['AR:NM_000044(UTR)'])
        self.assertEqual(define_transcripts(ar_utr),expected7)

        #Test when utr goes into exon 2
        hras_utr_exon2=self.exons['chr11'].search(int(534300), int(535450))
        expected8=(['HRAS'],['EXONIC'],['HRAS:NM_005343(UTR - exon 02)'])
        self.assertEqual(define_transcripts(hras_utr_exon2),expected8)

    def testGenomeIntervalTreeReverse(self):
        data=IntervalTree()
        for start,end,data in self.exons['chr17'].search(int(7577100)):
            if data['name']=='NM_000546':
                rev_exon=data
        for start,end,data in self.exons['chr17'].search(int(7577700)):
            if data['name']=='NM_000546':
                rev_intron=data

        #TP53:NM_000546, reverse strand, position 7577100 is in exon8,position 7577100 is in intron6
        self.assertEqual(rev_exon['exonNum'],'08')
        self.assertEqual(rev_intron['intronNum'],'06')

    def testGenomeIntervalTreeForward(self):
        data=IntervalTree()
        forward_exon=[]
        forward_intron=[]
        for start,end,data in self.exons['chrX'].search(int(66764988)):
            if data['name']=='NM_000044':
                forward_exon.append(data)
        for start,end,data in self.exons['chrX'].search(int(66942830)):
            if data['name']=='NM_000044':
                forward_intron.append(data)
        #AR:NM_000044, forward strand, position 66763880  is in exon1,position 66942830 is in intron7
        #Assert only 1 line of data is found
        self.assertEqual(len(forward_exon),1)
        self.assertEqual(len(forward_intron),1)
        self.assertEqual(forward_exon[0]['exonNum'],'01')
        self.assertEqual(forward_intron[0]['intronNum'],'07')

    def testRev5UTRstart(self):
        data=IntervalTree()
        rev_5_utr_start=[]
        for start,end,data in self.exons['chr17'].search(int(7571720)):
            if data['name']=='NM_000546':
                rev_5_utr_start.append(data.keys())
        self.assertEqual(len(rev_5_utr_start),1)        
        self.assertIn('UTR',rev_5_utr_start[0])

    def testRev5UTRstop(self):
        data=IntervalTree()
        rev_5_utr_stop=[]
        for start,end,data in self.exons['chr17'].search(int(7572925)):
            if data['name']=='NM_000546':
                rev_5_utr_stop.append(data.keys())
        self.assertEqual(len(rev_5_utr_stop),1)        
        self.assertIn('UTR',rev_5_utr_stop[0])

    def testRev3UTRstart(self):
        data=IntervalTree()
        ##3' UTR is  7579912-7590868
        rev_3_utr_start=[]
        for start,end,data in self.exons['chr17'].search(int(7579913)):
            if data['name']=='NM_000546':
                rev_3_utr_start.append(data.keys())
        self.assertEqual(len(rev_3_utr_start),1)
        self.assertIn('UTR',rev_3_utr_start[0])

    def testRev3UTRstart(self):
        data=IntervalTree()
        rev_3_utr_stop=[]
        for start,end,data in self.exons['chr17'].search(int(7590868)):
            if data['name']=='NM_000546':
                rev_3_utr_stop.append(data.keys())
        self.assertEqual(len(rev_3_utr_stop),1)
        self.assertIn('UTR',rev_3_utr_stop[0])
        
    def test3UTRstart(self):         
        data=IntervalTree()
        utr_3_start=[]
        #3' UTR is 66950461-66943683
        for start,end,data in self.exons['chrX'].search(int(66950461)):
            if data['name']=='NM_000044':
                utr_3_start.append(data.keys())
        self.assertEqual(len(utr_3_start),1)
        self.assertIn('UTR',utr_3_start[0])

    def test3UTRstop(self):
        data=IntervalTree()
        #Test 3' UTR
        utr_3_stop=[]
        for start,end,data in self.exons['chrX'].search(int(66943684)):
            if data['name']=='NM_000044':
                utr_3_stop.append(data.keys())
        self.assertEqual(len(utr_3_stop),1)
        self.assertIn('UTR',utr_3_stop[0])

    def test5UTRstart(self):        
        data=IntervalTree()
        utr_5_start=[]
        #5' UTR - 66763873-66764988
        for start,end,data in self.exons['chrX'].search(int(66763873)):
            if data['name']=='NM_000044':
                utr_5_start.append(data.keys())
        self.assertEqual(len(utr_5_start),1)
        self.assertIn('UTR',utr_5_start[0])

    def test5UTRstop(self):        
        data=IntervalTree()
        # Test 5' UTR
        utr_5_stop=[]
        for start,end,data in self.exons['chrX'].search(int(66764987)):
            if data['name']=='NM_000044':
                utr_5_stop.append(data.keys())
        self.assertEqual(len(utr_5_stop),1)
        self.assertIn('UTR',utr_5_stop[0])


    def testHRAS(self):
        data=IntervalTree()
        # Test UTR
        utr=[]
        for start,end,data in self.exons['chr11'].search(int(532242)):
            if data['name']=='NM_005343':
                utr.append(data.keys())
        self.assertEqual(len(utr),1)
        self.assertIn('UTR',utr[0])

        data=IntervalTree()
        # Test exon2
        exon2=[]
        for start,end,data in self.exons['chr11'].search(int(532750)):
            if data['name']=='NM_005343':
                exon2.append(data)
        self.assertEqual(len(exon2),1)
        self.assertEqual(exon2[0]['exonNum'],'05')

        # Test intron
        intron=[]
        for start,end,data in self.exons['chr11'].search(int(533000)):
            if data['name']=='NM_005343':
                intron.append(data)
        self.assertEqual(len(intron),1)
        self.assertEqual(intron[0]['intronNum'],'04')

        # Test UTR
        utr2=[]
        for start,end,data in self.exons['chr11'].search(int(532630)):
            if data['name']=='NM_005343':
                utr2.append(data.keys())
        self.assertEqual(len(utr2), 1)
        self.assertIn('UTR',utr2[0])
