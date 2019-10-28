'''
Test the annotsv_summary subcommand
'''
import subprocess
import filecmp
import StringIO
import sys
from os import path
import unittest
import logging
import pandas as pd
from munging.subcommands import annotsv_summary
from intervaltree import Interval, IntervalTree

from __init__ import TestBase
import __init__ as config
log = logging.getLogger(__name__)

annotsv_testfiles = path.join(config.datadir, 'annotsv')

class TestAnnotSV(TestBase):
    '''
    Test the annotsv_summary script
    '''
    
    def setUp(self):
        ''' Read test data into dataframe for use in various tests'''
        annotsv_df=pd.read_csv(path.join(annotsv_testfiles, 'small_annotsv.txt'), delimiter='\t', index_col=False, usecols=['SV chrom','SV start','SV end', 'ID', 'ALT','Gene name','NM','QUAL',
                                                                                                                            'FILTER','INFO','location','promoters','1000g_event', '1000g_max_AF', 
                                                                                                                            'Repeats_type_left', 'Repeats_type_right',
                                                                                                                            'DGV_GAIN_n_samples_with_SV','DGV_GAIN_n_samples_tested',
                                                                                                                            'DGV_LOSS_n_samples_with_SV','DGV_LOSS_n_samples_tested'])

        annotsv_df.fillna('', inplace=True)
        self.annotsv_df = annotsv_df

        bed_file = path.join(annotsv_testfiles, 'small_bed.bed')
        self.capture_trees = annotsv_summary.build_capture_trees(bed_file)

        flagged_fusions_file = path.join(annotsv_testfiles, 'small_flagged.txt')
        self.fusion_trees = annotsv_summary.parse_flagged_fusions_file(flagged_fusions_file)

    def testParseLength(self):
        '''Create leght if on same chrom,
        otherwiser print ITX
        '''
        annotsv_df=self.annotsv_df.copy()
        annotsv_df=annotsv_df.apply(annotsv_summary.parse_sv_event1, axis=1).apply(annotsv_summary.parse_sv_alt, axis=1).apply(annotsv_summary.parse_gene_promoter,axis=1).apply(annotsv_summary.parse_dgv, axis=1).apply(annotsv_summary.parse_repeats,axis=1).apply(annotsv_summary.parse_info, axis=1).apply(annotsv_summary.parse_location, axis=1)

        o_event='gridss137_4056o'
        input_o_data=self.annotsv_df.loc[(annotsv_df['ID']==o_event)]
        o_dict = annotsv_summary.collapse_event(input_o_data)
        o_event1 = annotsv_df.loc[annotsv_df['ID']==o_event,'Event1'].iloc[0]
        o_event2 = annotsv_df.loc[annotsv_df['ID']==o_event,'Event2'].iloc[0]
        event1=o_event1
        event2=o_event2
        length=annotsv_summary.parse_length(event1,event2)
        self.assertEqual(length, 1948852)

    def testParseLengthCTX(self):
        '''Create leght if on same chrom,
        otherwiser print CTX
        '''
        annotsv_df=self.annotsv_df.copy()
        annotsv_df=annotsv_df.apply(annotsv_summary.parse_sv_event1, axis=1).apply(annotsv_summary.parse_sv_alt, axis=1).apply(annotsv_summary.parse_gene_promoter,axis=1).apply(annotsv_summary.parse_dgv, axis=1).apply(annotsv_summary.parse_repeats,axis=1).apply(annotsv_summary.parse_info, axis=1).apply(annotsv_summary.parse_location, axis=1)

        o_event='gridss29_10153o'
        input_o_data=self.annotsv_df.loc[(annotsv_df['ID']==o_event)]
        o_dict = annotsv_summary.collapse_event(input_o_data)
        o_event1 = annotsv_df.loc[annotsv_df['ID']==o_event,'Event1'].iloc[0]
        o_event2 = annotsv_df.loc[annotsv_df['ID']==o_event,'Event2'].iloc[0]
        event1=o_event1
        event2=o_event2
        length=annotsv_summary.parse_length(event1,event2)
        self.assertEqual(length, 'CTX')

    def testParseQuality1(self):
        ''' Test quality parsing when all calls are above threshold
        '''
        annotsv_df=self.annotsv_df.copy()
        #filter all calls less than 200 quality
        annotsv_df=annotsv_summary.parse_quality(annotsv_df, quality=200)
        expected_quals=[366.85, 366.85, 366.85, 322.03, 322.03, 322.03, 322.03, 2268.1, 2268.1, 2268.1, 2268.1, 1002.51, 1002.51, 1002.51, 1002.51, 1184.9, 1184.9, 1184.9, 1184.9, 1109.12, 215.25, 215.25, 285.77, 285.77]
        quals=[x for x in annotsv_df['QUAL']]
        self.assertListEqual(sorted(quals), sorted(expected_quals))

    def testParseQuality2(self):
        ''' Test removing calls that do not meet quality threshold
        '''
        annotsv_df=self.annotsv_df.copy()
        #filter all calls less than 200 quality
        annotsv_df=annotsv_summary.parse_quality(annotsv_df, quality=250)
        expected_quals=[366.85, 366.85, 366.85, 322.03, 322.03, 322.03, 322.03, 2268.1, 2268.1, 2268.1, 2268.1, 1002.51, 1002.51, 1002.51, 1002.51, 1184.9, 1184.9, 1184.9, 1184.9, 1109.12]
        quals=[x for x in annotsv_df['QUAL']]
        self.assertListEqual(sorted(quals), sorted(expected_quals))

    def testParseSVALT(self):
        '''Parse the ALT breakend format into regular chr#:POS,
        returning Event2'''
        expected_input_alts=['A[7:55249011[', 'A[7:55249011[', ']GL:55248960]A', ']7:140490765]C', ']7:140490765]C', 'A[7:138541913[', 'A[7:138541913[', 'T[X:66766396[', 'T[X:66766396[', ']X:66766356]G', ']X:66766356]G', 'C[3:178921649[', 'C[3:178921649[', ']3:178921591]T', ']3:178921591]T', ']12:66451467]A', ']12:66451467]A', ']12:66451467]A', ']12:66451467]A', 'C[2:48028531[', 'T[7:98550704[', 'T[7:98550704[', ']7:98550669]T', ']7:98550669]T']
        #Make sure the input to the test hasn't changed
        df_alts=[x for x in self.annotsv_df['ALT']]
        self.assertListEqual(sorted(df_alts),sorted(expected_input_alts))

        #Make sure the function is working correctly
        event2_alt_df=self.annotsv_df.apply(annotsv_summary.parse_sv_alt, axis=1)
        #the 'if x==x' removes any 'nan' from this, which occurs for positions not in regular chrm1-23,X,Y
        event2_alts=[x for x in event2_alt_df['Event2'] if x==x]
        expected_event2_alts=['chr7:55249011', 'chr7:55249011', 'chr7:140490765', 'chr7:140490765', 'chr7:138541913', 'chr7:138541913', 'chrX:66766396', 'chrX:66766396', 'chrX:66766356', 'chrX:66766356', 'chr3:178921649', 'chr3:178921649', 'chr3:178921591', 'chr3:178921591', 'chr12:66451467', 'chr12:66451467', 'chr12:66451467', 'chr12:66451467', 'chr2:48028531', 'chr7:98550704', 'chr7:98550704', 'chr7:98550669', 'chr7:98550669']
        self.assertListEqual(sorted(event2_alts), sorted(expected_event2_alts))
        
    def testParseSVEvent1(self):
        ''' Combine fields to make Event1 '''
        event1_df=self.annotsv_df.apply(annotsv_summary.parse_sv_event1, axis=1)
        expected_event1=['chr7:55248960', 'chr7:55248960', 'chr7:55249011', 'chr7:138541913', 'chr7:138541913', 'chr7:140490765', 'chr7:140490765', 'chrX:66766356', 'chrX:66766356', 'chrX:66766396', 'chrX:66766396', 'chr3:178921591', 'chr3:178921591', 'chr3:178921649', 'chr3:178921649', 'chr2:48028531', 'chr2:48028531', 'chr2:48028531', 'chr2:48028531', 'chr12:66451467', 'chr7:98550671', 'chr7:98550671', 'chr7:98550704', 'chr7:98550704']
        event1=[x for x in event1_df['Event1']]
        self.assertListEqual(sorted(event1), sorted(expected_event1))

    def testParseInfo(self):
        '''Get the EventID for each read '''
        #Make sure the function is working correctly
        eventIDs_df=self.annotsv_df.apply(annotsv_summary.parse_info, axis=1)
        expected_eventIDs=['gridss129_14', 'gridss129_14', 'gridss129_14', 'gridss133_319', 'gridss133_319', 'gridss133_319', 'gridss133_319','gridss137_4056', 'gridss137_4056', 'gridss137_4056', 'gridss137_4056', 'gridss295_7', 'gridss295_7', 'gridss295_7', 'gridss295_7', 'gridss67_8', 'gridss67_8', 'gridss67_8', 'gridss67_8', 'gridss29_10153', 'gridss29_10153', 'gridss29_10153', 'gridss29_10153', 'gridss29_10153']
        eventIDs=[x for x in eventIDs_df['EventID']]
        self.assertListEqual(sorted(eventIDs), sorted(expected_eventIDs))

    def testParseGenePromoter(self):
         ''' Combine promoter and gene fields '''
         #Make sure the function is working correctly
         genes_df=self.annotsv_df.apply(annotsv_summary.parse_gene_promoter, axis=1)
         expected_genes=['EGFR/EGFR-AS1', 'EGFR', 'EGFR/EGFR-AS1', 'KIAA1549', 'KIAA1549', 'BRAF', 'BRAF', 'AR', 'AR', 'AR', 'AR', 'PIK3CA[Promoter]', 'PIK3CA[Promoter]', 'PIK3CA[Promoter]', 'PIK3CA[Promoter]', 'MSH6', 'MSH6', 'MSH6', 'MSH6', '', 'TRRAP', 'TRRAP', 'TRRAP', 'TRRAP']
         genes=[x for x in genes_df['Gene']]
         self.assertListEqual(sorted(genes), sorted(expected_genes))

    def testParseLocation(self):
         ''' Split location and remove duplicates'''

         #Make sure the function is working correctly
         locs_df=self.annotsv_df.apply(annotsv_summary.parse_location, axis=1)
         expected_locs=['intron18', 'intron16-intron6',  'intron8', 'exon1', 'exon1', 'intron5', 'intron5', 'intron3', 'intron3','intron37','intron37']
         locs=[x for x in locs_df['location'] if str(x) != 'nan' and str(x) !='']
         self.assertListEqual(sorted(locs), sorted(expected_locs))
    
    def testParseDGV(self):
        '''Combine DGV gain/loss columns into
        gain_n/gain_n_samples loss_n/loss_n_samples'''

        dgv_df=self.annotsv_df.apply(annotsv_summary.parse_dgv, axis=1)
        dgv_gain=[x for x in dgv_df['DGV_GAIN_found|tested']]
        dgv_lost=[x for x in dgv_df['DGV_LOSS_found|tested']]
        expected_dgv_gain=['0|0', '0|0', '0|0', '0|0','0|0', '0|0', '0|0', '0|0', '0|0', '0|0', '0|0', '0|0', '0|0', '0|0', '0|0', '11|33', '11|33', '11|33', '11|33', '0|0', '0|0', '0|0', '0|0', '0|0']
        expected_dgv_lost=['0|0', '0|0', '0|0', '0|0','0|0', '0|0', '0|0', '0|0', '0|0', '0|0', '0|0', '10|20', '10|20', '10|20', '10|20', '0|0', '0|0', '0|0', '0|0', '0|0', '0|0', '0|0', '0|0', '10|97']
        self.assertListEqual(sorted(dgv_gain),sorted(expected_dgv_gain))
        self.assertListEqual(sorted(dgv_lost),sorted(expected_dgv_lost))

    def testParseRepeats(self):
        ''' Combine left and right repeat info into one column'''
        #Make sure the function is working correctly
        repeats_df=self.annotsv_df.apply(annotsv_summary.parse_repeats, axis=1)
        expected_repeats=['MLT2B1[left];MLT2B1[right]','AluSx[left];AluSx[right]','(CGG)n[left];(CGG)n[right]','(CGG)n[left];(CGG)n[right]','AluJb[left];AluJb[right]','AluJb[left];AluJb[right]','(T)n/AluSx1[left];(T)n/AluSx1[right]','(TG)n[left];(TG)n[right]','MER4C/(TG)n[left];MER4C/(TG)n[right]']
        repeats=[x for x in repeats_df['Repeats'] if str(x) != 'nan' and str(x) !='']
        self.assertListEqual(sorted(repeats), sorted(expected_repeats))

    def testSmooshEventIntoOneLine(self):
        ''' Test Smooshing a multiline annotsv event into one line'''
        expected_result=['chr7:138541913','chr7:140490765','KIAA1549','BRAF','intron16-intron6','intron8-intron8', 1948852, 'NM_001354609;NM_001164665','322.03','LOW_QUAL','','','MLT2B1[left];MLT2B1[right]','AluSx[left];AluSx[right]','0|0','0|0']
        eventIDs_df=self.annotsv_df.apply(annotsv_summary.parse_sv_event1, axis=1).apply(annotsv_summary.parse_sv_alt, axis=1).apply(annotsv_summary.parse_gene_promoter,axis=1).apply(annotsv_summary.parse_dgv, axis=1).apply(annotsv_summary.parse_repeats,axis=1).apply(annotsv_summary.parse_info, axis=1)
        event='gridss137_4056'
        input_data=eventIDs_df.loc[(eventIDs_df['EventID']==event)]
        smooshed_result=annotsv_summary.smoosh_event_into_one_line(input_data.copy())
        self.assertListEqual(sorted(smooshed_result),sorted(expected_result))

    def testCollapseEvent(self):
        ''' Combine various column entries into one string each'''
        o_event='gridss137_4056o'
        input_o_data=self.annotsv_df.loc[(self.annotsv_df['ID']==o_event)]
        o_dict = annotsv_summary.collapse_event(input_o_data)
        expected_o_dict={'INFO': 'AS=1;ASQ=85.61;ASRP=3;ASSR=11;BA=0;BANRP=0;BANRPQ=0.00;BANSR=0;BANSRQ=0.00;BAQ=0.00;BASRP=0;BASSR=0;BEID=asm137-11351,asm137-6182;BEIDH=0,0;BEIDL=115,300;BQ=0.00;BSC=0;BSCQ=0.00;BUM=0;BUMQ=0.00;BVF=0;CAS=0;CASQ=0.00;CIPOS=-2,0;CIRPOS=-2,0;CQ=322.03;EVENT=gridss137_4056;HOMLEN=2;HOMSEQ=GA;IC=0;IHOMPOS=-2,0;IQ=0.00;PARID=gridss137_4056h;RAS=1;RASQ=140.29;REF=0;REFPAIR=0;RP=2;RPQ=21.04;SB=0.0;SC=1X1N1X176M;SR=4;SRQ=75.09;SVTYPE=BND;VF=9', 'Repeats_type_right': 'MLT2B1', 'Repeats_type_left': 'MLT2B1', 'NM': 'NM_001164665', 'DGV_GAIN_n_samples_tested': '0', 'DGV_LOSS_n_samples_with_SV': '0', 'promoters': '', '1000g_max_AF': '', 'DGV_GAIN_n_samples_with_SV': '0', 'SV end': '138541914', 'ID': 'gridss137_4056o', 'FILTER': 'LOW_QUAL', 'QUAL': '322.03', 'Gene name': 'KIAA1549', 'SV start': '138541913', 'DGV_LOSS_n_samples_tested': '0', 'ALT': ']7:140490765]C', '1000g_event': '', 'SV chrom': '7', 'location': 'intron16-intron6'}
        self.assertEqual(o_dict, expected_o_dict)

    def testFailures(self):
        ''' Test for the two types of failures: 
        1. missing o or h
        2. o_ref != h_alt
        '''
        testing_output=path.join(annotsv_testfiles, 'testing_output.txt')
        in_file=path.join(annotsv_testfiles, 'small_annotsv.txt')
        cmd=["./munge", "annotsv_summary",in_file, "-s", "-o", testing_output]
        failure=subprocess.check_output(cmd).split('\n')
        self.assertEqual('Calls did not match for events o gridss133_319o/h gridss133_319h, expected: o1 chr7:98550671 == h2 chr7:98550669; o2 chr7:98550704 == h1 chr7:98550704',failure[1])
        
    def testParseSingleton(self):
        annotsv_df=self.annotsv_df.copy()
        annotsv_df=annotsv_df.apply(annotsv_summary.parse_sv_event1, axis=1).apply(annotsv_summary.parse_sv_alt, axis=1).apply(annotsv_summary.parse_gene_promoter,axis=1).apply(annotsv_summary.parse_dgv, axis=1).apply(annotsv_summary.parse_repeats,axis=1).apply(annotsv_summary.parse_info, axis=1).apply(annotsv_summary.parse_location, axis=1)
        o_event='gridss133_319o'
        o_dict = annotsv_summary.collapse_event(annotsv_df.loc[(annotsv_df['ID']==o_event)])
        output=annotsv_summary.parse_singleton(o_dict)
        expected_output=['chr7:98550671', 'chr7:98550704', 'TRRAP', 'TRRAP', 'intron37', 'SINGLETON EVENT', 'NM_003496', '215.25', 'SINGLETON EVENT;LOW_QUAL', '', '', 'MER4C/(TG)n[left];MER4C/(TG)n[right]', 'SINGLETON EVENT', '0|0', '0|0']
        self.assertEqual(sorted(output), sorted(expected_output))

    def testParseBreakEnd(self):
        """Test parsing of a singleton event"""
        annotsv_df=self.annotsv_df.copy()
        annotsv_df=annotsv_df.apply(annotsv_summary.parse_sv_event1, axis=1).apply(annotsv_summary.parse_sv_alt, axis=1).apply(annotsv_summary.parse_gene_promoter,axis=1).apply(annotsv_summary.parse_dgv, axis=1).apply(annotsv_summary.parse_repeats,axis=1).apply(annotsv_summary.parse_info, axis=1).apply(annotsv_summary.parse_location, axis=1)
        o_event='gridss133_319o'
        o_dict = annotsv_summary.collapse_event(annotsv_df.loc[(annotsv_df['ID']==o_event)])
        output=annotsv_summary.parse_singleton(o_dict)
        expected_output=['chr7:98550671', 'chr7:98550704', 'TRRAP', 'TRRAP', 'intron37', 'SINGLETON EVENT', 'NM_003496', '215.25', 'SINGLETON EVENT;LOW_QUAL', '', '', 'MER4C/(TG)n[left];MER4C/(TG)n[right]', 'SINGLETON EVENT', '0|0', '0|0']
        self.assertEqual(sorted(output), sorted(expected_output))

    def testParseType1(self):
        """Test parsing of GENE_FUSION SVs"""
        expected_output=['FOO1','FOO2', 'chr1:1234567', 'GENE_FUSION']
        data = pd.DataFrame({'Gene1':['FOO1'], 'Gene2':['FOO2'], 'Event2': ['chr1:1234567']})
        data['Type'] = data.apply(annotsv_summary.parse_event_type, axis=1)
        output = data.iloc[0].tolist()
        self.assertEqual(sorted(output), sorted(expected_output))

    def testParseType2(self):
        """Test parsing of INTERGENIC SVs"""
        expected_output=['Intergenic','Intergenic','chr1:1234567', 'INTERGENIC']
        data = pd.DataFrame({'Gene1':['Intergenic'], 'Gene2':['Intergenic'], 'Event2': ['chr1:1234567']})
        data['Type'] = data.apply(annotsv_summary.parse_event_type, axis=1)
        output = data.iloc[0].tolist()
        self.assertEqual(sorted(output), sorted(expected_output))

    def testParseType3(self):
        """Test parsing of INTRAGENIC SVs with Promoter"""
        expected_output=['FOO1[Promoter]','FOO1', 'chr1:1234567', 'INTRAGENIC']
        data = pd.DataFrame({'Gene1':['FOO1[Promoter]'], 'Gene2':['FOO1'], 'Event2': ['chr1:1234567']})
        data['Type'] = data.apply(annotsv_summary.parse_event_type, axis=1)
        output = data.iloc[0].tolist()
        self.assertEqual(sorted(output), sorted(expected_output))
    
    def testParseType4(self):
        """Test parsing of INTRAGENIC with concatenated genes"""
        expected_output=['FOO1;FOO2','FOO2;FOO3', 'chr1:1234567', 'INTRAGENIC']
        data = pd.DataFrame({'Gene1':['FOO1;FOO2'], 'Gene2':['FOO2;FOO3'], 'Event2': ['chr1:1234567']})
        data['Type'] = data.apply(annotsv_summary.parse_event_type, axis=1)
        output = data.iloc[0].tolist()
        self.assertEqual(sorted(output), sorted(expected_output))

    def testParseType5(self):
        """Test parsing of SINGLETON SVs"""
        expected_output=['FOO1','FOO1', 'SingleBreakEnd', 'SINGLETON']
        data = pd.DataFrame({'Gene1':['FOO1'], 'Gene2':['FOO1'], 'Event2': ['SingleBreakEnd']})
        data['Type'] = data.apply(annotsv_summary.parse_event_type, axis=1)
        output = data.iloc[0].tolist()
        self.assertEqual(sorted(output), sorted(expected_output))

    def testBuildCaptureTree(self):
        """Test creation of collection of interval trees from BED file"""
        capture_trees = self.capture_trees.copy()
        self.assertEqual(len(capture_trees), 24)
        self.assertEqual(len(capture_trees['2']), 156)
        self.assertEqual(len(capture_trees['X']), 15)
        with self.assertRaises(KeyError):
          capture_trees['Z']
        self.assertEqual(capture_trees['2'][21226164].pop()[2], 'APOB')
        self.assertEqual(capture_trees['2'][21226200].pop()[2], 'APOB')
        self.assertEqual(capture_trees['2'][21226284].pop()[2], 'APOB')
        self.assertEqual(len(capture_trees['2'][24000000]), 0)

    def testParseCaptureIntent1(self):
        """Test capture intent parsing when Event1 is on-target and Event2 is off"""
        capture_trees = self.capture_trees.copy()
        expected_output=['chr2:21226164','chr2:10000000','YES']
        data = pd.DataFrame({'Event1':['chr2:21226164'], 'Event2':['chr2:10000000']})
        data['Intended_For_Capture'] = data.apply(annotsv_summary.parse_capture_intent, interval_trees=capture_trees, axis=1)
        output = data.iloc[0].tolist()
        self.assertEqual(sorted(output), sorted(expected_output))

    def testParseCaptureIntent2(self):
        """Test capture intent parsing when Event1 is off-target and Event2 is the last base in on-target interval"""
        capture_trees = self.capture_trees.copy()
        expected_output=['chr1:10000000','chr2:21226284','YES']
        data = pd.DataFrame({'Event1':['chr1:10000000'], 'Event2':['chr2:21226284']})
        data['Intended_For_Capture'] = data.apply(annotsv_summary.parse_capture_intent, interval_trees=capture_trees, axis=1)
        output = data.iloc[0].tolist()
        self.assertEqual(sorted(output), sorted(expected_output))
    
    def testParseCaptureIntent3(self):
        """Test capture intent parsing when both events are off-target"""
        capture_trees = self.capture_trees.copy()
        expected_output=['chr2:21226163','chr2:21226285','NO']
        data = pd.DataFrame({'Event1':['chr2:21226163'], 'Event2':['chr2:21226285']})
        data['Intended_For_Capture'] = data.apply(annotsv_summary.parse_capture_intent, interval_trees=capture_trees, axis=1)
        output = data.iloc[0].tolist()
        self.assertEqual(sorted(output), sorted(expected_output))

    def testParseFusionsFile(self):
        """Test creation of collection of interval trees from flagged fusions file"""
        fusion_trees = self.fusion_trees.copy()
        self.assertEqual(len(fusion_trees), 24)
        self.assertEqual(len(fusion_trees['chr2']), 2)
        self.assertEqual(len(fusion_trees['chrX']), 0)
        with self.assertRaises(KeyError):
          fusion_trees['chrZ']

        MSHA2 = fusion_trees['chr2'][47669646].pop()[2]
        self.assertEqual(len(MSHA2), 1)
        self.assertEqual(MSHA2['chr2'][38121356].pop()[2], 'BOLAND')
        with self.assertRaises(KeyError):
            MSHA2['chr7']

        BRAF = fusion_trees['chr7'][140433812].pop()[2]
        self.assertEqual(len(BRAF), 18)
        self.assertEqual(BRAF['chr1'][52556388].pop()[2], 'Quiver:BRAF-BTF3L4')
        with self.assertRaises(KeyError):
            BRAF['chr8']

    def testParseClinicalFusions1(self):
        """Test detection of Boland fusion"""
        fusion_trees = self.fusion_trees.copy()
        expected_output=['chr2:47669500', 'MSHA2', 'chr2:38121000', 'Intergenic', 'GENE_FUSION', 'BOLAND']
        data = pd.DataFrame({'Event1':['chr2:47669500'], 'Event2':['chr2:38121000'], 'Gene1':['MSHA2'], 'Gene2':['Intergenic']})
        data['Type'] = data.apply(annotsv_summary.parse_event_type, axis=1)
        data['Flagged_Fusions'] = data.apply(annotsv_summary.parse_clinical_fusions, fusion_partners=fusion_trees, axis=1)
        output = data.iloc[0].tolist()
        self.assertEqual(sorted(output), sorted(expected_output))

    def testParseClinicalFusions2(self):
        """Test detection of Quiver fusion"""
        fusion_trees = self.fusion_trees.copy()
        expected_output=['chr7:140433812', 'BRAF', 'chr3:11314009', 'ATG7', 'GENE_FUSION', '=HYPERLINK("http://quiver.archerdx.com/results?query=BRAF%3AATG7", "Quiver:BRAF-ATG7")']
        data = pd.DataFrame({'Event1':['chr7:140433812'], 'Event2':['chr3:11314009'], 'Gene1':['BRAF'], 'Gene2':['ATG7']})
        data['Type'] = data.apply(annotsv_summary.parse_event_type, axis=1)
        data['Flagged_Fusions'] = data.apply(annotsv_summary.parse_clinical_fusions, fusion_partners=fusion_trees, axis=1)
        output = data.iloc[0].tolist()
        self.assertEqual(sorted(output), sorted(expected_output))

    def testParseClinicalFusions3(self):
        """Test when not a flagged fusion"""
        fusion_trees = self.fusion_trees.copy()
        expected_output=['chr1:10000000', 'FOO1', 'chr2:10000000', 'FOO2', 'GENE_FUSION', 'NO']
        data = pd.DataFrame({'Event1':['chr1:10000000'], 'Event2':['chr2:10000000'], 'Gene1':['FOO1'], 'Gene2':['FOO2']})
        data['Type'] = data.apply(annotsv_summary.parse_event_type, axis=1)
        data['Flagged_Fusions'] = data.apply(annotsv_summary.parse_clinical_fusions, fusion_partners=fusion_trees, axis=1)
        output = data.iloc[0].tolist()
        self.assertEqual(sorted(output), sorted(expected_output))
