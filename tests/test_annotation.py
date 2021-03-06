"""
Test the annotation functions
"""

import logging
import os
import sys
import munging.annotation as ann
from intervaltree import IntervalTree
from __init__ import TestBase
import __init__ as config

log = logging.getLogger(__name__)

ann_testfiles = os.path.join(config.datadir, 'annotation')

class TestAnnotation(TestBase):

    def setUp(self):
        self.outdir = self.mkoutdir()
        # a small refgene table to test annotations
        self.refgene = os.path.join(ann_testfiles, 'ann_mini_refgene.tsv')
        # a GenomeIntervalTree made from that refgene
        self.gt=ann.GenomeIntervalTree.from_table(open(self.refgene, 'r'))
        # a collection of transcripts from a densely-packed section of chromosome 1
        self.transcripts=[x[2] for x in self.gt['1'][144146810:146467744]]

    def testGetLocation01(self):
        """
        Tests string arguments
        """
        loc = ann.get_location('1', '1', '100')
        self.assertEquals(loc, 'chr1:1-100')
        loc = ann.get_location('1', '1', '1')
        self.assertEquals(loc, 'chr1:1')

    def testGetLocation02(self):
        """
        Tests integer arguments
        """
        loc = ann.get_location(1, 1, 100)
        self.assertEquals(loc, 'chr1:1-100')
        loc = ann.get_location(1, 1, 1)
        self.assertEquals(loc, 'chr1:1')

    def testSplitString(self):
        """
        Tests spliting a string given string of split characters
        Used to split path names in crawlers
        """
        result=ann.multi_split('/home/genetics/data/run_info', '/_')
        self.assertEquals(result, ['home','genetics','data','run','info'])


    def testSplitChrLoc(self):
        """
        Tests spliting the chr_loc from chr1:1-100 to 1, 1, 100
        Returns strings, not integers
        """
        result01=ann.split_chr_loc('chr1:1-100 ')
        result02=ann.split_chr_loc('1:1-100')
        self.assertEquals(result01, ('1','1','100'))
        self.assertEquals(result02, ('1','1','100'))

    def testSplitStringinTwo(self):
        """
        Tests spliting a string into two on | or ,
        Return -1, -1 if input is None
        Used to split columns of data in summary
        """
        result01=ann.split_string_in_two("1|99")
        result02=ann.split_string_in_two(None)
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
        result01,ref_reads01,varreads01=ann.build_variant_id(data1)
        result02,ref_reads02,varreads02=ann.build_variant_id(data2)
        self.assertListEqual([result01,ref_reads01,varreads01],['X_12321_12321_A_T','10','11'])
        self.assertListEqual([result02,ref_reads02,varreads02], ['X_1234_1256_G_C','12','13'])

    def testFixPfx(self):
        self.assertEqual(ann.fix_pfx('48_A03_BROv7-HA0186-NA12878'), '48_A03_BROv7_HA0186_NA12878')
        self.assertEqual(ann.fix_pfx('48_A03_BROv7-HA0186-NA12878 '), '48_A03_BROv7_HA0186_NA12878')
        self.assertEqual(ann.fix_pfx('UNK-124-455'), 'UNK_124_455')
        self.assertEqual(ann.fix_pfx('LMG-240'), 'LMG240')

    def testTranscriptHRAS(self):
        """
        HRAS is a gene on the - strand of chromosome 11; the NM_005343 transcript has
        six exons, of which Exon 1 and Exon 6 are non-coding.
        """
        # some facts about HRAS
        prefix = 'HRAS:NM_005343'
        chrom = '11'
        length = 3326
        num_exons = 6
        num_coding_exons = 4
        tx_start = 532241
        tx_end = 535567 - 1
        cd_start = 532635
        cd_end = 534322 - 1
        first_coding_exon_start = 532630
        first_coding_exon_end = 532755 - 1
        last_coding_exon_start = 534211
        last_coding_exon_end = 534375 - 1
        first_intron_start = 532522
        first_intron_end = 532630 - 1
        last_intron_start = 534375
        last_intron_end = 535415 - 1
        
        # get the HRAS transcript
        t = self.gt[chrom][tx_start:tx_end+1].pop()[2]

        # test its properties
        self.assertEqual(str(t), prefix)
        self.assertEqual(len(t), length)

        ### test the get_annotation() function ###
        # annotating the whole transcript
        self.assertEqual(t.get_annotation(tx_start, tx_end + 1), prefix + '(UTR - UTR)')
        self.assertEqual(t.get_annotation(tx_start, tx_end + 1, report_utr=False), prefix + '(exon 01 - exon 06)')
        self.assertEqual(t.get_annotation(0, sys.maxint), prefix + '(UTR - UTR)')
        # annotating outside the transcript
        self.assertEqual(t.get_annotation(0), None)
        self.assertEqual(t.get_annotation(0, tx_start), None)
        self.assertEqual(t.get_annotation(tx_end + 1, sys.maxint), None)
        # the 5' UTR
        self.assertEqual(t.get_annotation(tx_start), prefix + '(UTR)')
        self.assertEqual(t.get_annotation(tx_start, report_utr=False), prefix + '(exon 06)')
        # the 3' UTR
        self.assertEqual(t.get_annotation(tx_end), prefix + '(UTR)')
        self.assertEqual(t.get_annotation(tx_end, report_utr=False), prefix + '(exon 01)')
        # the boundary between UTR and the start of coding
        self.assertEqual(t.get_annotation(cd_start - 1), prefix + '(UTR)')
        self.assertEqual(t.get_annotation(cd_start, report_utr=False), prefix + '(exon 05)')
        self.assertEqual(t.get_annotation(cd_start), prefix + '(exon 05)')
        self.assertEqual(t.get_annotation(cd_start, report_utr=False), prefix + '(exon 05)')
        # the boundary between the end of coding and the UTR
        self.assertEqual(t.get_annotation(cd_end + 1), prefix + '(UTR)')
        self.assertEqual(t.get_annotation(cd_end + 1, report_utr=False), prefix + '(exon 02)')
        self.assertEqual(t.get_annotation(cd_end), prefix + '(exon 02)')
        self.assertEqual(t.get_annotation(cd_end, report_utr=False), prefix + '(exon 02)')
        # the entire first coding exon
        self.assertEqual(t.get_annotation(first_coding_exon_start, first_coding_exon_end + 1), prefix + '(exon 05 - UTR)')
        self.assertEqual(t.get_annotation(first_coding_exon_start, first_coding_exon_end + 1, report_utr=False), prefix + '(exon 05)')
        # the entire last coding exon
        self.assertEqual(t.get_annotation(last_coding_exon_start, last_coding_exon_end + 1), prefix + '(UTR - exon 02)')
        self.assertEqual(t.get_annotation(last_coding_exon_start, last_coding_exon_end + 1, report_utr=False), prefix + '(exon 02)')
        # tx_start to cd_end
        self.assertEqual(t.get_annotation(tx_start, cd_end + 1), prefix + '(exon 02 - UTR)')
        self.assertEqual(t.get_annotation(tx_start, cd_end + 1, report_utr=False), prefix + '(exon 02 - exon 06)')
        # cd_start to tx_end
        self.assertEqual(t.get_annotation(cd_start, tx_end + 1), prefix + '(UTR - exon 05)')
        self.assertEqual(t.get_annotation(cd_start, tx_end + 1, report_utr=False), prefix + '(exon 01 - exon 05)')
        # an intron between non-coding portions of exons
        self.assertEqual(t.get_annotation(first_intron_start, first_intron_end + 1), prefix + '(intron 05)')
        self.assertEqual(t.get_annotation(first_intron_start, first_intron_end + 1, report_utr=False), prefix + '(intron 05)')
        self.assertEqual(t.get_annotation(last_intron_start, last_intron_end + 1), prefix + '(intron 01)')
        self.assertEqual(t.get_annotation(last_intron_start, last_intron_end + 1, report_utr=False), prefix + '(intron 01)')

        ### test the get_exons() function ###
        # all the exons
        self.assertEqual(len(t.get_exons()), num_coding_exons)
        self.assertEqual(len(t.get_exons(report_utr=False)), num_exons)
        self.assertEqual(len(t.get_exons(tx_start, tx_end)), num_coding_exons)
        self.assertEqual(len(t.get_exons(tx_start, tx_end, report_utr=False)), num_exons)
        # no exons
        self.assertEqual(len(t.get_exons(0, tx_start, report_utr=False)), 0)
        self.assertEqual(len(t.get_exons(tx_end + 1, sys.maxint, report_utr=False)), 0)
        self.assertEqual(len(t.get_exons(0, cd_start, report_utr=True)), 0)
        self.assertEqual(len(t.get_exons(cd_end + 1, sys.maxint, report_utr=True)), 0)
        # barely one exon
        self.assertEqual(len(t.get_exons(0, tx_start + 1, report_utr=False)), 1)
        self.assertEqual(len(t.get_exons(0, cd_start + 1, report_utr=True)), 1)
        self.assertEqual(len(t.get_exons(tx_end, sys.maxint, report_utr=False)), 1)
        self.assertEqual(len(t.get_exons(cd_end, sys.maxint, report_utr=True)), 1)

        ## test the get region_types() function ###
        # whole transcript
        self.assertEqual(t.get_region_types(tx_start, tx_end + 1), {'EXONIC', 'INTRONIC', 'UTR'})
        self.assertEqual(t.get_region_types(tx_start, tx_end + 1, report_utr=False), {'EXONIC', 'INTRONIC'})
        self.assertEqual(t.get_region_types(0, sys.maxint), {'EXONIC', 'INTRONIC', 'UTR'})
        # outside the transcript
        self.assertEqual(t.get_region_types(0), set())
        self.assertEqual(t.get_region_types(sys.maxint), set())
        self.assertEqual(t.get_region_types(0, tx_start), set())
        self.assertEqual(t.get_region_types(tx_end + 1, sys.maxint), set())
        # the 5' UTR
        self.assertEqual(t.get_region_types(tx_start), {'UTR'})
        self.assertEqual(t.get_region_types(tx_start, report_utr=False), {'EXONIC'})
        # the 3' UTR
        self.assertEqual(t.get_region_types(tx_end), {'UTR'})
        self.assertEqual(t.get_region_types(tx_end, report_utr=False), {'EXONIC'})
        # the boundary between coding and non-coding at the 5' end
        self.assertEqual(t.get_region_types(cd_start - 1), {'UTR'})
        self.assertEqual(t.get_region_types(cd_start, report_utr=False), {'EXONIC'})
        self.assertEqual(t.get_region_types(cd_start), {'EXONIC'})
        self.assertEqual(t.get_region_types(cd_start, report_utr=False), {'EXONIC'})
        # the boundary between coding and non-coding at the 3' end
        self.assertEqual(t.get_region_types(cd_end + 1), {'UTR'})
        self.assertEqual(t.get_region_types(cd_end + 1, report_utr=False), {'EXONIC'})
        self.assertEqual(t.get_region_types(cd_end), {'EXONIC'})
        self.assertEqual(t.get_region_types(cd_end, report_utr=False), {'EXONIC'})
        # an intron between non-coding portions of exons
        self.assertEqual(t.get_region_types(first_intron_start, first_intron_end + 1), {'INTRONIC'})
        self.assertEqual(t.get_region_types(first_intron_start, first_intron_end + 1, report_utr=False), {'INTRONIC'})
        self.assertEqual(t.get_region_types(last_intron_start, last_intron_end + 1), {'INTRONIC'})
        self.assertEqual(t.get_region_types(last_intron_start, last_intron_end + 1, report_utr=False), {'INTRONIC'})


    def testTranscriptCHEK1(self):
        """
        CHEK1 is a gene on the + strand of chromosome 11; the NM_001114121 transcript has
        14 exons, of which Exon 1 and Exon 14 are non-coding.
        """
        # some facts about CHEK1
        prefix = 'CHEK1:NM_001114121'
        chrom = '11'
        length = 51120
        num_exons = 14
        num_coding_exons = 12
        tx_start = 125495030
        tx_end = 125546150 - 1
        cd_start = 125496663
        cd_end = 125525215 - 1
        first_coding_exon_start = 125496643
        first_coding_exon_end = 125496728 - 1
        last_coding_exon_start = 125525119
        last_coding_exon_end = 125525242 - 1
        first_intron_start = 125495907
        first_intron_end = 125496643 - 1
        last_intron_start = 125525242
        last_intron_end = 125545822 - 1
        
        # get the CHEK1 transcript
        t = self.gt[chrom][tx_start:tx_end+1].pop()[2]

        # test its properties
        self.assertEqual(str(t), prefix)
        self.assertEqual(len(t), length)

        ### test the get_annotation() function ###
        # annotating the whole transcript
        self.assertEqual(t.get_annotation(tx_start, tx_end + 1), prefix + '(UTR - UTR)')
        self.assertEqual(t.get_annotation(tx_start, tx_end + 1, report_utr=False), prefix + '(exon 01 - exon 14)')
        self.assertEqual(t.get_annotation(0, sys.maxint), prefix + '(UTR - UTR)')
        # annotating outside the transcript
        self.assertEqual(t.get_annotation(0), None)
        self.assertEqual(t.get_annotation(0, tx_start), None)
        self.assertEqual(t.get_annotation(tx_end + 1, sys.maxint), None)
        # the 5' UTR
        self.assertEqual(t.get_annotation(tx_start), prefix + '(UTR)')
        self.assertEqual(t.get_annotation(tx_start, report_utr=False), prefix + '(exon 01)')
        # the 3' UTR
        self.assertEqual(t.get_annotation(tx_end), prefix + '(UTR)')
        self.assertEqual(t.get_annotation(tx_end, report_utr=False), prefix + '(exon 14)')
        # the boundary between UTR and the start of coding
        self.assertEqual(t.get_annotation(cd_start - 1), prefix + '(UTR)')
        self.assertEqual(t.get_annotation(cd_start, report_utr=False), prefix + '(exon 02)')
        self.assertEqual(t.get_annotation(cd_start), prefix + '(exon 02)')
        self.assertEqual(t.get_annotation(cd_start, report_utr=False), prefix + '(exon 02)')
        # the boundary between the end of coding and the UTR
        self.assertEqual(t.get_annotation(cd_end + 1), prefix + '(UTR)')
        self.assertEqual(t.get_annotation(cd_end + 1, report_utr=False), prefix + '(exon 13)')
        self.assertEqual(t.get_annotation(cd_end), prefix + '(exon 13)')
        self.assertEqual(t.get_annotation(cd_end, report_utr=False), prefix + '(exon 13)')
        # the entire first coding exon
        self.assertEqual(t.get_annotation(first_coding_exon_start, first_coding_exon_end + 1), prefix + '(UTR - exon 02)')
        self.assertEqual(t.get_annotation(first_coding_exon_start, first_coding_exon_end + 1, report_utr=False), prefix + '(exon 02)')
        # the entire last coding exon
        self.assertEqual(t.get_annotation(last_coding_exon_start, last_coding_exon_end + 1), prefix + '(exon 13 - UTR)')
        self.assertEqual(t.get_annotation(last_coding_exon_start, last_coding_exon_end + 1, report_utr=False), prefix + '(exon 13)')
        # tx_start to cd_end
        self.assertEqual(t.get_annotation(tx_start, cd_end + 1), prefix + '(UTR - exon 13)')
        self.assertEqual(t.get_annotation(tx_start, cd_end + 1, report_utr=False), prefix + '(exon 01 - exon 13)')
        # cd_start to tx_end
        self.assertEqual(t.get_annotation(cd_start, tx_end + 1), prefix + '(exon 02 - UTR)')
        self.assertEqual(t.get_annotation(cd_start, tx_end + 1, report_utr=False), prefix + '(exon 02 - exon 14)')
        # an intron between non-coding portions of exons
        self.assertEqual(t.get_annotation(first_intron_start, first_intron_end + 1), prefix + '(intron 01)')
        self.assertEqual(t.get_annotation(first_intron_start, first_intron_end + 1, report_utr=False), prefix + '(intron 01)')
        self.assertEqual(t.get_annotation(last_intron_start, last_intron_end + 1), prefix + '(intron 13)')
        self.assertEqual(t.get_annotation(last_intron_start, last_intron_end + 1, report_utr=False), prefix + '(intron 13)')

        ### test the get_exons() function ###
        # all the exons
        self.assertEqual(len(t.get_exons()), num_coding_exons)
        self.assertEqual(len(t.get_exons(report_utr=False)), num_exons)
        self.assertEqual(len(t.get_exons(tx_start, tx_end)), num_coding_exons)
        self.assertEqual(len(t.get_exons(tx_start, tx_end, report_utr=False)), num_exons)
        # no exons
        self.assertEqual(len(t.get_exons(0, tx_start, report_utr=False)), 0)
        self.assertEqual(len(t.get_exons(tx_end + 1, sys.maxint, report_utr=False)), 0)
        self.assertEqual(len(t.get_exons(0, cd_start, report_utr=True)), 0)
        self.assertEqual(len(t.get_exons(cd_end + 1, sys.maxint, report_utr=True)), 0)
        # barely one exon
        self.assertEqual(len(t.get_exons(0, tx_start + 1, report_utr=False)), 1)
        self.assertEqual(len(t.get_exons(0, cd_start + 1, report_utr=True)), 1)
        self.assertEqual(len(t.get_exons(tx_end, sys.maxint, report_utr=False)), 1)
        self.assertEqual(len(t.get_exons(cd_end, sys.maxint, report_utr=True)), 1)

        ## test the get region_types() function ###
        # whole transcript
        self.assertEqual(t.get_region_types(tx_start, tx_end + 1), {'EXONIC', 'INTRONIC', 'UTR'})
        self.assertEqual(t.get_region_types(tx_start, tx_end + 1, report_utr=False), {'EXONIC', 'INTRONIC'})
        self.assertEqual(t.get_region_types(0, sys.maxint), {'EXONIC', 'INTRONIC', 'UTR'})
        # outside the transcript
        self.assertEqual(t.get_region_types(0), set())
        self.assertEqual(t.get_region_types(sys.maxint), set())
        self.assertEqual(t.get_region_types(0, tx_start), set())
        self.assertEqual(t.get_region_types(tx_end + 1, sys.maxint), set())
        # the 5' UTR
        self.assertEqual(t.get_region_types(tx_start), {'UTR'})
        self.assertEqual(t.get_region_types(tx_start, report_utr=False), {'EXONIC'})
        # the 3' UTR
        self.assertEqual(t.get_region_types(tx_end), {'UTR'})
        self.assertEqual(t.get_region_types(tx_end, report_utr=False), {'EXONIC'})
        # the boundary between coding and non-coding at the 5' end
        self.assertEqual(t.get_region_types(cd_start - 1), {'UTR'})
        self.assertEqual(t.get_region_types(cd_start, report_utr=False), {'EXONIC'})
        self.assertEqual(t.get_region_types(cd_start), {'EXONIC'})
        self.assertEqual(t.get_region_types(cd_start, report_utr=False), {'EXONIC'})
        # the boundary between coding and non-coding at the 3' end
        self.assertEqual(t.get_region_types(cd_end + 1), {'UTR'})
        self.assertEqual(t.get_region_types(cd_end + 1, report_utr=False), {'EXONIC'})
        self.assertEqual(t.get_region_types(cd_end), {'EXONIC'})
        self.assertEqual(t.get_region_types(cd_end, report_utr=False), {'EXONIC'})
        # an intron between non-coding portions of exons
        self.assertEqual(t.get_region_types(first_intron_start, first_intron_end + 1), {'INTRONIC'})
        self.assertEqual(t.get_region_types(first_intron_start, first_intron_end + 1, report_utr=False), {'INTRONIC'})
        self.assertEqual(t.get_region_types(last_intron_start, last_intron_end + 1), {'INTRONIC'})
        self.assertEqual(t.get_region_types(last_intron_start, last_intron_end + 1, report_utr=False), {'INTRONIC'})

    def testTranscriptTACSTD2(self):
        """
        TACSTD2 is a gene on the - strand of chromosome 1; the NM_002353 transcript has one exon.
        """
        # some facts about TACSTD2
        prefix = 'TACSTD2:NM_002353'
        chrom = '1'
        length = 2072
        num_exons = 1
        num_coding_exons = 1
        tx_start = 59041094
        tx_end = 59043166 - 1
        cd_start = 59041856
        cd_end = 59042828 - 1

        # get the CHEK1 transcript
        t = self.gt[chrom][tx_start:tx_end+1].pop()[2]

        # test its properties
        self.assertEqual(str(t), prefix)
        self.assertEqual(len(t), length)

        ### test the get_annotation() function ###
        # annotating the whole transcript
        self.assertEqual(t.get_annotation(tx_start, tx_end + 1), prefix + '(UTR - UTR)')
        self.assertEqual(t.get_annotation(tx_start, tx_end + 1, report_utr=False), prefix + '(exon 01)')
        self.assertEqual(t.get_annotation(0, sys.maxint), prefix + '(UTR - UTR)')
        # annotating outside the transcript
        self.assertEqual(t.get_annotation(0), None)
        self.assertEqual(t.get_annotation(0, tx_start), None)
        self.assertEqual(t.get_annotation(tx_end + 1, sys.maxint), None)
        # the 5' UTR
        self.assertEqual(t.get_annotation(tx_start), prefix + '(UTR)')
        self.assertEqual(t.get_annotation(tx_start, report_utr=False), prefix + '(exon 01)')
        # the 3' UTR
        self.assertEqual(t.get_annotation(tx_end), prefix + '(UTR)')
        self.assertEqual(t.get_annotation(tx_end, report_utr=False), prefix + '(exon 01)')
        # the boundary between UTR and the start of coding
        self.assertEqual(t.get_annotation(cd_start - 1), prefix + '(UTR)')
        self.assertEqual(t.get_annotation(cd_start, report_utr=False), prefix + '(exon 01)')
        self.assertEqual(t.get_annotation(cd_start), prefix + '(exon 01)')
        self.assertEqual(t.get_annotation(cd_start, report_utr=False), prefix + '(exon 01)')
        # the boundary between the end of coding and the UTR
        self.assertEqual(t.get_annotation(cd_end + 1), prefix + '(UTR)')
        self.assertEqual(t.get_annotation(cd_end + 1, report_utr=False), prefix + '(exon 01)')
        self.assertEqual(t.get_annotation(cd_end), prefix + '(exon 01)')
        self.assertEqual(t.get_annotation(cd_end, report_utr=False), prefix + '(exon 01)')
        # tx_start to cd_end
        self.assertEqual(t.get_annotation(tx_start, cd_end + 1), prefix + '(exon 01 - UTR)')
        self.assertEqual(t.get_annotation(tx_start, cd_end + 1, report_utr=False), prefix + '(exon 01)')
        # cd_start to tx_end
        self.assertEqual(t.get_annotation(cd_start, tx_end + 1), prefix + '(UTR - exon 01)')
        self.assertEqual(t.get_annotation(cd_start, tx_end + 1, report_utr=False), prefix + '(exon 01)')

        ### test the get_exons() function ###
        # all the exons
        self.assertEqual(len(t.get_exons()), num_coding_exons)
        self.assertEqual(len(t.get_exons(report_utr=False)), num_exons)
        self.assertEqual(len(t.get_exons(tx_start, tx_end)), num_coding_exons)
        self.assertEqual(len(t.get_exons(tx_start, tx_end, report_utr=False)), num_exons)
        # no exons
        self.assertEqual(len(t.get_exons(0, tx_start, report_utr=False)), 0)
        self.assertEqual(len(t.get_exons(tx_end + 1, sys.maxint, report_utr=False)), 0)
        self.assertEqual(len(t.get_exons(0, cd_start, report_utr=True)), 0)
        self.assertEqual(len(t.get_exons(cd_end + 1, sys.maxint, report_utr=True)), 0)
        # barely one exon
        self.assertEqual(len(t.get_exons(0, tx_start + 1, report_utr=False)), 1)
        self.assertEqual(len(t.get_exons(0, cd_start + 1, report_utr=True)), 1)
        self.assertEqual(len(t.get_exons(tx_end, sys.maxint, report_utr=False)), 1)
        self.assertEqual(len(t.get_exons(cd_end, sys.maxint, report_utr=True)), 1)

        ## test the get region_types() function ###
        # whole transcript
        self.assertEqual(t.get_region_types(tx_start, tx_end + 1), {'EXONIC', 'UTR'})
        self.assertEqual(t.get_region_types(tx_start, tx_end + 1, report_utr=False), {'EXONIC'})
        self.assertEqual(t.get_region_types(0, sys.maxint), {'EXONIC', 'UTR'})
        # outside the transcript
        self.assertEqual(t.get_region_types(0), set())
        self.assertEqual(t.get_region_types(sys.maxint), set())
        self.assertEqual(t.get_region_types(0, tx_start), set())
        self.assertEqual(t.get_region_types(tx_end + 1, sys.maxint), set())
        # the 5' UTR
        self.assertEqual(t.get_region_types(tx_start), {'UTR'})
        self.assertEqual(t.get_region_types(tx_start, report_utr=False), {'EXONIC'})
        # the 3' UTR
        self.assertEqual(t.get_region_types(tx_end), {'UTR'})
        self.assertEqual(t.get_region_types(tx_end, report_utr=False), {'EXONIC'})
        # the boundary between coding and non-coding at the 5' end
        self.assertEqual(t.get_region_types(cd_start - 1), {'UTR'})
        self.assertEqual(t.get_region_types(cd_start, report_utr=False), {'EXONIC'})
        self.assertEqual(t.get_region_types(cd_start), {'EXONIC'})
        self.assertEqual(t.get_region_types(cd_start, report_utr=False), {'EXONIC'})
        # the boundary between coding and non-coding at the 3' end
        self.assertEqual(t.get_region_types(cd_end + 1), {'UTR'})
        self.assertEqual(t.get_region_types(cd_end + 1, report_utr=False), {'EXONIC'})
        self.assertEqual(t.get_region_types(cd_end), {'EXONIC'})
        self.assertEqual(t.get_region_types(cd_end, report_utr=False), {'EXONIC'})

    def testGeneInfoFromTranscripts(self):
        """Test extracting gene names from overlapping transcripts on both strands"""
        all_genes = ['LIX1L',
                    'LOC100996724',
                    'NBPF10',
                    'NBPF20',
                    'NBPF25P',
                    'NBPF8',
                    'NBPF9',
                    'NOTCH2NL',
                    'PDE4DIP',
                    'PDZK1',
                    'RNF115',
                    'SEC22B']
        # all the transcripts
        self.assertEqual(ann.gene_info_from_transcripts(self.transcripts), all_genes)
        self.assertEqual(ann.gene_info_from_transcripts(self.transcripts, 0, sys.maxint), all_genes)
        # none of the transcripts
        self.assertEqual(ann.gene_info_from_transcripts(self.transcripts, 0), ['Intergenic'])
        self.assertEqual(ann.gene_info_from_transcripts(self.transcripts, 0, 144146810), ['Intergenic'])
        self.assertEqual(ann.gene_info_from_transcripts(self.transcripts, 146467744 + 1, sys.maxint), ['Intergenic'])
        # overlapping transcripts on both strands
        self.assertEqual(ann.gene_info_from_transcripts(self.transcripts, 145076000), ['NBPF20', 'NBPF9', 'PDE4DIP'])

    def testRegionInfoFromTranscripts(self):
        """Test extracting transcript region types from overlapping transcripts on both strands"""
        # all the transcripts
        self.assertEqual(ann.region_info_from_transcripts(self.transcripts, 0, sys.maxint), ['EXONIC', 'INTRONIC', 'UTR']) 
        self.assertEqual(ann.region_info_from_transcripts(self.transcripts, 0, sys.maxint, report_utr=False), ['EXONIC', 'INTRONIC']) 
        # none of the transcripts
        self.assertEqual(ann.region_info_from_transcripts(self.transcripts, 0), ['Intergenic'])
        self.assertEqual(ann.region_info_from_transcripts(self.transcripts, 0, 144146810), ['Intergenic'])
        self.assertEqual(ann.region_info_from_transcripts(self.transcripts, 146467744 + 1, sys.maxint), ['Intergenic'])
        # overlapping transcripts on both strands
        self.assertEqual(ann.region_info_from_transcripts(self.transcripts, 145076000), ['INTRONIC', 'UTR'])
        self.assertEqual(ann.region_info_from_transcripts(self.transcripts, 145076000, report_utr=False), ['EXONIC', 'INTRONIC'])

    def testTranscriptInfoFromTranscripts(self):
        """Test extracting transcript annotations from overlapping transcripts on both strands"""
        all_annotations = ['LIX1L:NM_153713(UTR - UTR)',
                           'LOC100996724:NR_144516(UTR - UTR)',
                           'LOC100996724:NR_144517(UTR - UTR)',
                           'NBPF10:NM_001039703(UTR - exon 86)',
                           'NBPF10:NM_001302371(UTR - exon 90)',
                           'NBPF20:NM_001278267(UTR - UTR)',
                           'NBPF25P:NR_104217(UTR - UTR)',
                           'NBPF8:NM_001037501(UTR - UTR)',
                           'NBPF8:NR_102404(UTR - UTR)',
                           'NBPF8:NR_102405(UTR - UTR)',
                           'NBPF9:NM_001037675(UTR - UTR)',
                           'NBPF9:NM_001277444(UTR - UTR)',
                           'NOTCH2NL:NM_203458(UTR - UTR)',
                           'PDE4DIP:NM_001002810(UTR - UTR)',
                           'PDE4DIP:NM_001002812(UTR - UTR)',
                           'PDE4DIP:NM_001195260(UTR - UTR)',
                           'PDE4DIP:NM_001195261(UTR - UTR)',
                           'PDE4DIP:NM_001198832(UTR - UTR)',
                           'PDE4DIP:NM_001198834(UTR - UTR)',
                           'PDE4DIP:NM_014644(UTR - UTR)',
                           'PDE4DIP:NM_022359(UTR - UTR)',
                           'PDZK1:NM_001201325(UTR - UTR)',
                           'PDZK1:NM_001201326(UTR - UTR)',
                           'PDZK1:NM_002614(UTR - UTR)',
                           'RNF115:NM_014455(UTR - UTR)',
                           'SEC22B:NM_004892(UTR - UTR)']
        all_annotations_no_utr = ['LIX1L:NM_153713(exon 01 - exon 06)',
                                  'LOC100996724:NR_144516(exon 01 - exon 45)',
                                  'LOC100996724:NR_144517(exon 01 - exon 43)',
                                  'NBPF10:NM_001039703(exon 01 - exon 86)',
                                  'NBPF10:NM_001302371(exon 01 - exon 90)',
                                  'NBPF20:NM_001278267(exon 01 - exon 131)',
                                  'NBPF25P:NR_104217(exon 01 - exon 18)',
                                  'NBPF8:NM_001037501(exon 01 - exon 21)',
                                  'NBPF8:NR_102404(exon 01 - exon 19)',
                                  'NBPF8:NR_102405(exon 01 - exon 18)',
                                  'NBPF9:NM_001037675(exon 01 - exon 22)',
                                  'NBPF9:NM_001277444(exon 01 - exon 26)',
                                  'NOTCH2NL:NM_203458(exon 01 - exon 05)',
                                  'PDE4DIP:NM_001002810(exon 01 - exon 04)',
                                  'PDE4DIP:NM_001002812(exon 01 - exon 23)',
                                  'PDE4DIP:NM_001195260(exon 01 - exon 05)',
                                  'PDE4DIP:NM_001195261(exon 01 - exon 05)',
                                  'PDE4DIP:NM_001198832(exon 01 - exon 46)', 
                                  'PDE4DIP:NM_001198834(exon 01 - exon 44)',
                                  'PDE4DIP:NM_014644(exon 01 - exon 44)',
                                  'PDE4DIP:NM_022359(exon 01 - exon 07)',
                                  'PDZK1:NM_001201325(exon 01 - exon 09)',
                                  'PDZK1:NM_001201326(exon 01 - exon 07)',
                                  'PDZK1:NM_002614(exon 01 - exon 10)',
                                  'RNF115:NM_014455(exon 01 - exon 09)',
                                  'SEC22B:NM_004892(exon 01 - exon 06)']
        # all the transcripts
        self.assertEqual(ann.transcript_info_from_transcripts(self.transcripts, 0, sys.maxint), all_annotations) 
        self.assertEqual(ann.transcript_info_from_transcripts(self.transcripts, 0, sys.maxint, report_utr=False), all_annotations_no_utr) 
        # none of the transcripts
        self.assertEqual(ann.transcript_info_from_transcripts(self.transcripts, 0), [])
        self.assertEqual(ann.transcript_info_from_transcripts(self.transcripts, 0, 144146810), [])
        self.assertEqual(ann.transcript_info_from_transcripts(self.transcripts, 146467744 + 1, sys.maxint), [])
        # overlapping transcripts on both strands
        self.assertEqual(ann.transcript_info_from_transcripts(self.transcripts, 145076000),
                         ['NBPF20:NM_001278267(intron 53)', 'NBPF9:NM_001277444(intron 21)', 'PDE4DIP:NM_022359(UTR)'])
        self.assertEqual(ann.transcript_info_from_transcripts(self.transcripts, 145076000, report_utr=False),
                         ['NBPF20:NM_001278267(intron 53)', 'NBPF9:NM_001277444(intron 21)', 'PDE4DIP:NM_022359(exon 01)'])
