"""
Test the filter functions
"""

import os
from os import path
import unittest
import logging
import pprint
import sys
import json
import csv

from operator import itemgetter
from itertools import ifilter
from collections import namedtuple, defaultdict
from munging import filters 
from munging.utils import Path, walker
from munging import parsers


from __init__ import TestBase
import __init__ as config

log = logging.getLogger(__name__)

testfiles = path.join(config.datadir, '101010_HA0000_OncoPlex1','output')
qualityfiles = path.join(config.datadir, 'quality_metrics')

class TestParsers(TestBase):
    """Test each of the parsers are returning the correct fieldnames, 
    prefixes and variant_key list"""
    def setUp(self):
        self.outdir = self.mkoutdir()

    def tearDown(self):
        pass

    def testSNPParser(self):
        """
        Test for correct fieldname parsing 
        """
        specimens = defaultdict(dict)
        annotation = {} 
        prefixes = []
        variant_keys = []
        sort_order=['5437_E05_OPXv4_NA12878_MA0013','6037_E05_OPXv4_NA12878_HA0201','0228T_CON_OPXv4_INT']
        analysis_type='parsers.parse_snp'
        files = ifilter(filters.any_analysis, walker(testfiles))  
        chosen_parser='{}(files, specimens, annotation, prefixes, variant_keys, sort_order)'.format(analysis_type)
        specimens, annotation, prefixes, fieldnames, variant_keys=eval(chosen_parser)
        self.assertListEqual(prefixes,['5437_NA12878_Ref|Var','6037_NA12878_Ref|Var','0228T_Ref|Var', 'Count'])
        self.assertListEqual(fieldnames, ['Position', 'Ref_Base', 'Var_Base', 'Gene', 'Variant_Type',
                                          'Transcripts', 'Clinically_Flagged', 'Cosmic', 'Segdup', 
                                          'Polyphen', 'Sift', 'Mutation_Taster', 'Gerp', 
                                          'UW_Freq', 'UW_Count', 'UW_DEC_p',
                                          '1000g_ALL', 'EVS_esp6500_ALL', '1000g_AMR', 'EVS_esp6500_AA', '1000g_EUR',
                                          'EVS_esp6500_EU', '1000g_SAS','1000g_EAS', '1000g_AFR', 'ADA_Alter_Splice', 'RF_Alter_Splice', 
                                          '5437_NA12878_Ref|Var','6037_NA12878_Ref|Var','0228T_Ref|Var', 'Count'])
        self.assertListEqual(variant_keys, ['Position', 'Ref_Base', 'Var_Base'])
        
    def testCNVGeneParser(self):
        specimens = defaultdict(dict)
        annotation = {} 
        prefixes = []
        variant_keys = []
        sort_order=['0228T_CON_OPXv4_INT','5437_E05_OPXv4_NA12878_MA0013','6037_E05_OPXv4_NA12878_HA0201']
        files = ifilter(filters.any_analysis, walker(testfiles))  
        analysis_type='parsers.parse_cnv_gene'
        chosen_parser='{}(files, specimens, annotation, prefixes, variant_keys, sort_order)'.format(analysis_type)
        specimens, annotation, prefixes, fieldnames, variant_keys=eval(chosen_parser)
        self.assertListEqual(prefixes,['0228T_Log', '5437_NA12878_Log', '6037_NA12878_Log'])
        self.assertListEqual(fieldnames, ['Position', 'Gene', 'Transcripts', '0228T_Log', '5437_NA12878_Log', '6037_NA12878_Log'])
        self.assertListEqual(variant_keys, ['Position', 'Gene'])

    def testCNVExonParser(self):
        specimens = defaultdict(dict)
        annotation = {} 
        prefixes = []
        variant_keys = []
        sort_order=['0228T_CON_OPXv4_INT','5437_E05_OPXv4_NA12878_MA0013','6037_E05_OPXv4_NA12878_HA0201']
        analysis_type='parsers.parse_cnv_exon'
        files = ifilter(filters.any_analysis, walker(testfiles))  
        chosen_parser='{}(files, specimens, annotation, prefixes, variant_keys, sort_order)'.format(analysis_type)
        specimens, annotation, prefixes, fieldnames, variant_keys=eval(chosen_parser)
        self.assertListEqual(prefixes,['0228T_Log', '5437_NA12878_Log', '6037_NA12878_Log'])
        self.assertListEqual(fieldnames, ['Position', 'Gene', 'Transcripts', '0228T_Log', '5437_NA12878_Log', '6037_NA12878_Log'])
        self.assertListEqual(variant_keys, ['Position', 'Gene'])
    
    def testQualityParser(self):
        specimens = defaultdict(dict)
        annotation = {} 
        prefixes = []
        variant_keys = []
        sort_order=['0228T_CON_OPXv4_INT','5437_E05_OPXv4_NA12878_MA0013','6037_E05_OPXv4_NA12878_HA0201']
        analysis_type='parsers.parse_quality'
        files = ifilter(filters.any_analysis, walker(testfiles))  
        chosen_parser='{}(files, specimens, annotation, prefixes, variant_keys, sort_order)'.format(analysis_type)
        specimens, annotation, prefixes, fieldnames, variant_keys=eval(chosen_parser)
        self.assertListEqual(fieldnames, ['MEAN_TARGET_COVERAGE', '0228T','5437_NA12878','6037_NA12878'])
        self.assertListEqual(variant_keys, ['MEAN_TARGET_COVERAGE'])
        self.assertListEqual(prefixes,['0228T','5437_NA12878','6037_NA12878'])
       
    def testPindelParser(self):
        specimens = defaultdict(dict)
        annotation = {} 
        prefixes = []
        variant_keys = []
        sort_order=['0228T_CON_OPXv4_INT','5437_E05_OPXv4_NA12878_MA0013','6037_E05_OPXv4_NA12878_HA0201']
        analysis_type='parsers.parse_pindel'
        files = ifilter(filters.any_analysis, walker(testfiles))  
        chosen_parser='{}(files, specimens, annotation, prefixes, variant_keys, sort_order)'.format(analysis_type)
        specimens, annotation, prefixes, fieldnames, variant_keys=eval(chosen_parser)
        self.assertListEqual(prefixes,['0228T', '5437_NA12878', '6037_NA12878','Count'])
        self.assertListEqual(fieldnames, ['Position', 'Gene', 'Size', 'Gene_Region', 'Event_Type','Transcripts', '0228T', '5437_NA12878', '6037_NA12878','Count'])
        self.assertListEqual(variant_keys, ['Position', 'Gene', 'Size'])
        
    def testClinFlaggedParser(self):
        specimens = defaultdict(dict)
        annotation = {} 
        prefixes = []
        variant_keys = []
        sort_order=['0228T_CON_OPXv4_INT','5437_E05_OPXv4_NA12878_MA0013','6037_E05_OPXv4_NA12878_HA0201']
        analysis_type='parsers.parse_clin_flagged'
        files = ifilter(filters.any_analysis, walker(testfiles))  
        chosen_parser='{}(files, specimens, annotation, prefixes, variant_keys, sort_order)'.format(analysis_type)
        specimens, annotation, prefixes, fieldnames, variant_keys=eval(chosen_parser)
        self.assertListEqual(prefixes,['0228T_Variants', '5437_NA12878_Variants', '6037_NA12878_Variants'])
        self.assertListEqual(fieldnames, ['Position', 'Ref_Base', 'Var_Base', 'Clinically_Flagged', '0228T_Variants', '5437_NA12878_Variants', '6037_NA12878_Variants'])
        self.assertListEqual(variant_keys, ['Position', 'Ref_Base', 'Var_Base'])
        
    def testHSParser(self):
        variant_keys = []
        files = ifilter(filters.hs_file_finder, walker(qualityfiles))
        fname=open(path.join(qualityfiles, '6037_E05_OPXv4_NA12878_HA0201.hs_metrics'),'rU')
        lines=fname.readlines()
        output_dict, variant_keys = parsers.parse_hsmetrics(lines, variant_keys)
        self.assertListEqual(sorted(variant_keys),sorted(['PF_READS',
                                                          'PF_UNIQUE_READS',
                                                          'PCT_PF_UQ_READS',
                                                          'PF_UQ_READS_ALIGNED',
                                                          'PCT_SELECTED_BASES',
                                                          'PCT_OFF_BAIT',
                                                          'MEAN_TARGET_COVERAGE',
                                                          'PCT_USABLE_BASES_ON_TARGET',
                                                          'ZERO_CVG_TARGETS_PCT',
                                                          'AT_DROPOUT',
                                                          'GC_DROPOUT']))
        self.assertDictContainsSubset({'MEAN_TARGET_COVERAGE':'614.820203'}, output_dict)

    def testQualityMetricsParser(self):
        variant_keys = []
        files = ifilter(filters.quality_file_finder, walker(qualityfiles))
        fname=open(path.join(qualityfiles, '6037_E05_OPXv4_NA12878_HA0201.quality_metrics'),'rU')
        lines=fname.readlines()
        output_dict, variant_keys = parsers.parse_qualitymetrics(lines, variant_keys)
        self.assertListEqual(sorted(variant_keys),sorted(['UNPAIRED_READS_EXAMINED',
                                                          'READ_PAIRS_EXAMINED',
                                                          'UNMAPPED_READS',  
                                                          'UNPAIRED_READ_DUPLICATES',
                                                          'READ_PAIR_DUPLICATES',
                                                          'READ_PAIR_OPTICAL_DUPLICATES',    
                                                          'PERCENT_DUPLICATION',     
                                                          'ESTIMATED_LIBRARY_SIZE']))
        self.assertDictContainsSubset({'PERCENT_DUPLICATION':'0.130625'}, output_dict)

    def testMSIFlaggedParser(self):
        specimens = defaultdict(dict)
        annotation = {} 
        prefixes = []
        variant_keys = []
        sort_order=['6037_E05_OPXv4_NA12878_HA0201','0228T_CON_OPXv4_INT','5437_E05_OPXv4_NA12878_MA0013']
        analysis_type='parsers.parse_msi_flagged'
        files = ifilter(filters.any_analysis, walker(testfiles))  
        chosen_parser='{}(files, specimens, annotation, prefixes, variant_keys, sort_order)'.format(analysis_type)
        specimens, annotation, prefixes, fieldnames, variant_keys=eval(chosen_parser)
        self.assertListEqual(prefixes,['6037_NA12878_Variants|Total', '6037_NA12878_Status', '0228T_Variants|Total', '0228T_Status', '5437_NA12878_Variants|Total', '5437_NA12878_Status'])
        self.assertListEqual(fieldnames, ['Position', 'Ref_Base', 'Var_Base', 'Clinically_Flagged', '6037_NA12878_Variants|Total', '6037_NA12878_Status', '0228T_Variants|Total', '0228T_Status', '5437_NA12878_Variants|Total', '5437_NA12878_Status'])
        self.assertListEqual(variant_keys, ['Position', 'Ref_Base', 'Var_Base'])
        self.assertEqual(specimens[('chr7:55259524','T','A')]['0228T_Status'], 'REVIEW') #less than 100 reads
        self.assertEqual(specimens[('chr3:37034946', 'G', 'A')]['0228T_Status'], 'POS')
        self.assertEqual(specimens[('chr12:25380283', 'C', 'T')]['0228T_Status'], 'NEG')
        self.assertEqual(specimens[('chr13:32936674', 'C', 'T')]['0228T_Status'], 'IND')

    def testGLTFlaggedParser(self):
        specimens = defaultdict(dict)
        annotation = {} 
        prefixes = []
        variant_keys = []
        sort_order=['6037_E05_OPXv4_NA12878_HA0201','0228T_CON_OPXv4_INT','5437_E05_OPXv4_NA12878_MA0013']
        analysis_type='parsers.parse_glt_flagged'
        files = ifilter(filters.any_analysis, walker(testfiles))  
        chosen_parser='{}(files, specimens, annotation, prefixes, variant_keys,sort_order)'.format(analysis_type)
        specimens, annotation, prefixes, fieldnames, variant_keys=eval(chosen_parser)
        self.assertListEqual(prefixes,['6037_NA12878_Variants|Total', '6037_NA12878_Status', '0228T_Variants|Total', '0228T_Status', '5437_NA12878_Variants|Total', '5437_NA12878_Status'])
        self.assertListEqual(fieldnames, ['Position', 'Ref_Base', 'Var_Base', 'Clinically_Flagged', '6037_NA12878_Variants|Total', '6037_NA12878_Status', '0228T_Variants|Total', '0228T_Status', '5437_NA12878_Variants|Total', '5437_NA12878_Status'])
        self.assertListEqual(variant_keys, ['Position', 'Ref_Base', 'Var_Base'])
        self.assertEqual(specimens[('chr7:55259524','T','A')]['0228T_Status'], 'REVIEW')#less than 100 reads
        self.assertEqual(specimens[('chr3:37034946', 'G', 'A')]['0228T_Status'], 'HET') #.4-.65 
        self.assertEqual(specimens[('chr2:215661788','C','T')]['0228T_Status'], 'HOMO')#>=.98
        self.assertEqual(specimens[('chr13:32936674', 'C', 'T')]['0228T_Status'], 'NEG')#<.1

    def testHotSpotFlaggedParser(self):
        specimens = defaultdict(dict)
        annotation = {} 
        prefixes = []
        variant_keys = []
        sort_order=['6037_E05_OPXv4_NA12878_HA0201','0228T_CON_OPXv4_INT','5437_E05_OPXv4_NA12878_MA0013']
        analysis_type='parsers.parse_hotspot_flagged'
        files = ifilter(filters.any_analysis, walker(testfiles))  
        chosen_parser='{}(files, specimens, annotation, prefixes, variant_keys,sort_order)'.format(analysis_type)
        specimens, annotation, prefixes, fieldnames, variant_keys=eval(chosen_parser)
        self.assertListEqual(prefixes,['6037_NA12878_Variants|Total', '6037_NA12878_Status', '0228T_Variants|Total', '0228T_Status', '5437_NA12878_Variants|Total', '5437_NA12878_Status'])
        self.assertListEqual(fieldnames, ['Position', 'Ref_Base', 'Var_Base', 'Clinically_Flagged', '6037_NA12878_Variants|Total', '6037_NA12878_Status', '0228T_Variants|Total', '0228T_Status', '5437_NA12878_Variants|Total', '5437_NA12878_Status'])
        self.assertListEqual(variant_keys, ['Position', 'Ref_Base', 'Var_Base'])
        self.assertEqual(specimens[('chr7:55259524','T','A')]['0228T_Status'], 'REVIEW')#less than 100 reads
        self.assertEqual(specimens[('chr3:37034946', 'G', 'A')]['0228T_Status'], 'HET')#.2-.7
        self.assertEqual(specimens[('chr2:215661788','C','T')]['0228T_Status'], 'HOMO')#>.7
        self.assertEqual(specimens[('chr13:32936674', 'C', 'T')]['0228T_Status'], 'NEG')#<.1

    def testAnnotSVParser(self):
        specimens = defaultdict(dict)
        annotation = {} 
        prefixes = []
        variant_keys = []
        sort_order=['0228T_CON_OPXv4_INT','5437_E05_OPXv4_NA12878_MA0013','6037_E05_OPXv4_NA12878_HA0201']
        analysis_type='parsers.parse_annotsv'
        files = ifilter(filters.any_analysis, walker(testfiles))  
        chosen_parser='{}(files, specimens, annotation, prefixes, variant_keys, sort_order)'.format(analysis_type)
        specimens, annotation, prefixes, fieldnames, variant_keys=eval(chosen_parser)
        self.assertListEqual(prefixes,['0228T', '5437_NA12878', '6037_NA12878','Count'])
        self.assertListEqual(fieldnames, ['Event1','Event2','Gene1','Gene2','NM','0228T', '5437_NA12878', '6037_NA12878','Count'])
        self.assertListEqual(variant_keys, ['Event1','Event2','Gene1','Gene2','NM'])
        self.assertEqual(len(specimens), 19)

    def testAmpliconParser(self):
        specimens = defaultdict(dict)
        annotation = {} 
        prefixes = []
        variant_keys = []
        sort_order=['0228T_CON_OPXv4_INT','5437_E05_OPXv4_NA12878_MA0013','6037_E05_OPXv4_NA12878_HA0201']
        analysis_type='parsers.parse_amplicon'
        files = ifilter(filters.any_analysis, walker(testfiles))
        chosen_parser='{}(files, specimens, annotation, prefixes, variant_keys, sort_order)'.format(analysis_type)
        specimens, annotation, prefixes, fieldnames, variant_keys=eval(chosen_parser)
        self.assertListEqual(prefixes,['0228T', '5437_NA12878', '6037_NA12878'])
        self.assertListEqual(fieldnames, ['Position', 'Probe','0228T', '5437_NA12878', '6037_NA12878'])
        self.assertListEqual(variant_keys, ['Position', 'Probe'])
        #Should heave 753 entries for hhv3, the only assay running this parser
        self.assertEqual(len(specimens), 753)
