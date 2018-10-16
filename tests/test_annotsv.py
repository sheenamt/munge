'''
Test the annotsv_summary subcommand
'''
from os import path
import unittest
import logging
import pandas as pd
from munging.subcommands import annotsv_summary

from __init__ import TestBase
import __init__ as config
log = logging.getLogger(__name__)

annotsv_testfiles = path.join(config.datadir, 'annotsv', 'small_annotsv.txt')

class TestAnnotSV(TestBase):
    '''
    Test the annotsv_summary script
    '''
    
    def setUp(self):
        ''' Read test data into dataframe for use in various tests'''
        annotsv_df=pd.read_csv(annotsv_testfiles, delimiter='\t', index_col=False, usecols=['SV chrom','SV start','SV end', 'ID', 'ALT','Gene name','NM','QUAL',
                                                                                             'FILTER','INFO','location','promoters','1000g_event', '1000g_max_AF', 
                                                                                             'Repeats_type_left', 'Repeats_type_right',
                                                                                             'DGV_GAIN_n_samples_with_SV','DGV_GAIN_n_samples_tested',
                                                                                             'DGV_LOSS_n_samples_with_SV','DGV_LOSS_n_samples_tested'])
        annotsv_df.fillna('', inplace=True)
        self.annotsv_df = annotsv_df

    def testParseSVALT(self):
        '''Parse the ALT breakend format into regular chr#:POS,
        returning Event2'''

        #Make sure the input to the test hasn't changed
        expected_input_alts=['A[7:55249011[', 'A[7:55249011[', ']7:55248960]A', ']7:140490765]C', ']7:140490765]C', 'A[7:138541913[', 'A[7:138541913[', 'T[X:66766396[', 'T[X:66766396[', ']X:66766356]G', ']X:66766356]G', 'C[3:178921649[', 'C[3:178921649[', ']3:178921591]T', ']3:178921591]T', ']12:66451467]A', ']12:66451467]A', ']12:66451467]A', ']12:66451467]A', 'C[2:48028531[']
        df_alts=[x for x in self.annotsv_df['ALT']]
        self.assertListEqual(sorted(df_alts),sorted(expected_input_alts))

        #Make sure the function is working correctly
        event2_alt_df=self.annotsv_df.apply(annotsv_summary.parse_sv_alt, axis=1)
        expected_event2_alts=['chr7:55249011', 'chr7:55249011', 'chr7:55248960', 'chr7:140490765', 'chr7:140490765', 'chr7:138541913', 'chr7:138541913', 'chrX:66766396', 'chrX:66766396', 'chrX:66766356', 'chrX:66766356', 'chr3:178921649', 'chr3:178921649', 'chr3:178921591', 'chr3:178921591', 'chr12:66451467', 'chr12:66451467', 'chr12:66451467', 'chr12:66451467', 'chr2:48028531']
        event2_alts=[x for x in event2_alt_df['Event2']]
        self.assertListEqual(sorted(event2_alts), sorted(expected_event2_alts))
        
    def testParseSVEvent1(self):
        ''' Combine fields to make Event1 
        '''

        event1_df=self.annotsv_df.apply(annotsv_summary.parse_sv_event1, axis=1)
        expected_event1=['chr7:55248960', 'chr7:55248960', 'chr7:55249011', 'chr7:138541913', 'chr7:138541913', 'chr7:140490765', 'chr7:140490765', 'chrX:66766356', 'chrX:66766356', 'chrX:66766396', 'chrX:66766396', 'chr3:178921591', 'chr3:178921591', 'chr3:178921649', 'chr3:178921649', 'chr2:48028531', 'chr2:48028531', 'chr2:48028531', 'chr2:48028531', 'chr12:66451467']
        event1=[x for x in event1_df['Event1']]
        self.assertListEqual(sorted(event1), sorted(expected_event1))


    def testParseInfo(self):
        '''Get the EventID for each read '''

        #Make sure the function is working correctly
        eventIDs_df=self.annotsv_df.apply(annotsv_summary.parse_info, axis=1)
        expected_eventIDs=['gridss129_14', 'gridss129_14', 'gridss129_14', 'gridss137_4056', 'gridss137_4056', 'gridss137_4056', 'gridss137_4056', 'gridss295_7', 'gridss295_7', 'gridss295_7', 'gridss295_7', 'gridss67_8', 'gridss67_8', 'gridss67_8', 'gridss67_8', 'gridss29_10153', 'gridss29_10153', 'gridss29_10153', 'gridss29_10153', 'gridss29_10153']
        eventIDs=[x for x in eventIDs_df['EventID']]
        self.assertListEqual(sorted(eventIDs), sorted(expected_eventIDs))

    def testParseGenePromoter(self):
         ''' Combine promoter and gene fields '''

         #Make sure the function is working correctly
         genes_df=self.annotsv_df.apply(annotsv_summary.parse_gene_promoter, axis=1)
         expected_genes=['EGFR/EGFR-AS1', 'EGFR', 'EGFR/EGFR-AS1', 'KIAA1549', 'KIAA1549', 'BRAF', 'BRAF', 'AR', 'AR', 'AR', 'AR', 'PIK3CA[Promoter]', 'PIK3CA[Promoter]', 'PIK3CA[Promoter]', 'PIK3CA[Promoter]', 'MSH6', 'MSH6', 'MSH6', 'MSH6', '']
         genes=[x for x in genes_df['Gene']]
         self.assertListEqual(sorted(genes), sorted(expected_genes))

    def testParseLocation(self):
         ''' Split location and remove duplicates'''

         #Make sure the function is working correctly
         locs_df=self.annotsv_df.apply(annotsv_summary.parse_location, axis=1)
         expected_locs=['intron18', 'intron16-intron6',  'intron8', 'exon1', 'exon1', 'intron5', 'intron5', 'intron3', 'intron3']
         locs=[x for x in locs_df['location'] if str(x) != 'nan' and str(x) !='']
         self.assertListEqual(sorted(locs), sorted(expected_locs))
    
    def testParseDGV(self):
        '''Combine DGV gain/loss columns into
        gain_n/gain_n_samples loss_n/loss_n_samples'''

        dgv_df=self.annotsv_df.apply(annotsv_summary.parse_dgv, axis=1)
        dgv_gain=[x for x in dgv_df['DGV_GAIN_found|tested']]
        dgv_lost=[x for x in dgv_df['DGV_LOSS_found|tested']]
        expected_dgv_gain=['0|0', '0|0', '0|0', '0|0', '0|0', '0|0', '0|0', '0|0', '0|0', '0|0', '0|0', '11|33', '11|33', '11|33', '11|33', '0|0', '0|0', '0|0', '0|0', '0|0']
        expected_dgv_lost=['0|0', '0|0', '0|0', '0|0', '0|0', '0|0', '0|0', '10|20', '10|20', '10|20', '10|20', '0|0', '0|0', '0|0', '0|0', '0|0', '0|0', '0|0', '0|0', '10|97']
        self.assertListEqual(sorted(dgv_gain),sorted(expected_dgv_gain))
        self.assertListEqual(sorted(dgv_lost),sorted(expected_dgv_lost))

    def testParseRepeats(self):
        ''' Combine left and right repeat info into one column'''

        #Make sure the function is working correctly
        repeats_df=self.annotsv_df.apply(annotsv_summary.parse_repeats, axis=1)
        expected_repeats=['MLT2B1[left];MLT2B1[right]','AluSx[left];AluSx[right]','(CGG)n[left];(CGG)n[right]','(CGG)n[left];(CGG)n[right]','AluJb[left];AluJb[right]','AluJb[left];AluJb[right]','(T)n/AluSx1[left];(T)n/AluSx1[right]']
        repeats=[x for x in repeats_df['Repeats'] if str(x) != 'nan' and str(x) !='']
        self.assertListEqual(sorted(repeats), sorted(expected_repeats))

    def testSmooshEventIntoOneLine(self):
        ''' Test Smooshing a multiline annotsv event into one line'''

        expected_result=['chr7:138541913','chr7:140490765','KIAA1549','BRAF','intron16-intron6','intron8-intron8','NM_001354609;NM_001164665','322.03','LOW_QUAL','','','MLT2B1[left];MLT2B1[right]','AluSx[left];AluSx[right]','0|0','0|0']
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


