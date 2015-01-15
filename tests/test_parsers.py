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

testfiles = path.join(config.datadir, 'analysis_files')
testMSIfile = path.join(config.datadir, 'MSI')

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
        files = ifilter(filters.any_analysis, walker(testfiles))  
        analysis_type='parsers.parse_snp'
        chosen_parser='{}(files, specimens, annotation, prefixes, variant_keys)'.format(analysis_type)
        specimens, annotation, prefixes, fieldnames, variant_keys=eval(chosen_parser)   
        self.assertListEqual(sorted(prefixes),sorted(['5437_NA12878_Ref|Var','6037_NA12878_Ref|Var','0228T_Ref|Var', 'Count']))
        self.assertListEqual(sorted(fieldnames), sorted(['Position', 'Ref_Base', 'Var_Base', 'Gene', 'Variant_Type',
                                         'Transcripts', 'Clinically_Flagged', 'Cosmic', 'Segdup', 
                                         'Polyphen', 'Sift', 'Mutation_Taster', 'Gerp', 'HiSeq_Freq',
                                         'HiSeq_Count', 'MiSeq_Freq', 'MiSeq_Count', '1000g_ALL', 
                                         'EVS_esp6500_ALL', '1000g_AMR', 'EVS_esp6500_AA', '1000g_EUR',
                                         'EVS_esp6500_EU', '1000g_ASN', '1000g_AFR', '6037_NA12878_Ref|Var',
                                                         '5437_NA12878_Ref|Var', '0228T_Ref|Var', 'Count']))
        self.assertListEqual(variant_keys, ['Position', 'Ref_Base', 'Var_Base'])
        
    def testCNVGeneParser(self):
        specimens = defaultdict(dict)
        annotation = {} 
        prefixes = []
        variant_keys = []
        files = ifilter(filters.any_analysis, walker(testfiles))  
        analysis_type='parsers.parse_cnv_gene'
        chosen_parser='{}(files, specimens, annotation, prefixes, variant_keys)'.format(analysis_type)
        specimens, annotation, prefixes, fieldnames, variant_keys=eval(chosen_parser)   
        self.assertListEqual(sorted(prefixes),sorted(['0228T_Log', '5437_NA12878_Log', '6037_NA12878_Log']))
        self.assertListEqual(sorted(fieldnames), sorted(['Position', 'Gene', 'Transcripts', '0228T_Log', '5437_NA12878_Log', '6037_NA12878_Log']))
        self.assertListEqual(variant_keys, ['Position', 'Gene'])

    def testCNVExonParser(self):
        specimens = defaultdict(dict)
        annotation = {} 
        prefixes = []
        variant_keys = []
        files = ifilter(filters.any_analysis, walker(testfiles))  
        analysis_type='parsers.parse_cnv_exon'
        chosen_parser='{}(files, specimens, annotation, prefixes, variant_keys)'.format(analysis_type)
        specimens, annotation, prefixes, fieldnames, variant_keys=eval(chosen_parser)   
        self.assertListEqual(sorted(prefixes),sorted(['0228T_Log', '5437_NA12878_Log', '6037_NA12878_Log']))
        self.assertListEqual(sorted(fieldnames), sorted(['Position', 'Gene', 'Transcripts', '0228T_Log', '5437_NA12878_Log', '6037_NA12878_Log']))
        self.assertListEqual(variant_keys, ['Position', 'Gene'])
    
    def testQualityParser(self):
        specimens = defaultdict(dict)
        annotation = {} 
        prefixes = []
        variant_keys = []
        files = ifilter(filters.any_analysis, walker(testfiles))
        analysis_type='parsers.parse_quality'
        chosen_parser='{}(files, specimens, annotation, prefixes, variant_keys)'.format(analysis_type)
        specimens, annotation, prefixes, fieldnames, variant_keys=eval(chosen_parser)        
        self.assertListEqual(fieldnames, ['MEAN_TARGET_COVERAGE', '0228T','5437_NA12878','6037_NA12878'])
        self.assertListEqual(variant_keys, ['MEAN_TARGET_COVERAGE'])
        self.assertListEqual(sorted(prefixes),sorted(['0228T','5437_NA12878','6037_NA12878']))
       
    def testPindelParser(self):
        specimens = defaultdict(dict)
        annotation = {} 
        prefixes = []
        variant_keys = []
        files = ifilter(filters.any_analysis, walker(testfiles))  
        analysis_type='parsers.parse_pindel'
        chosen_parser='{}(files, specimens, annotation, prefixes, variant_keys)'.format(analysis_type)
        specimens, annotation, prefixes, fieldnames, variant_keys=eval(chosen_parser)      
        self.assertListEqual(sorted(prefixes),sorted(['0228T', '5437_NA12878', '6037_NA12878','Count']))
        self.assertListEqual(sorted(fieldnames), sorted(['Position', 'Gene', 'Gene_Region', 'Event_Type', 'Size', 'Transcripts', '0228T', '5437_NA12878', '6037_NA12878','Count']))
        self.assertListEqual(variant_keys, ['Position', 'Gene'])
        
    def testClinFlaggedParser(self):
        specimens = defaultdict(dict)
        annotation = {} 
        prefixes = []
        variant_keys = []
        files = ifilter(filters.any_analysis, walker(testfiles))  
        analysis_type='parsers.parse_clin_flagged'
        chosen_parser='{}(files, specimens, annotation, prefixes, variant_keys)'.format(analysis_type)
        specimens, annotation, prefixes, fieldnames, variant_keys=eval(chosen_parser)    
        self.assertListEqual(sorted(prefixes),sorted(['0228T_Reads', '5437_NA12878_Reads', '6037_NA12878_Reads']))
        self.assertListEqual(sorted(fieldnames), sorted(['Position', 'Ref_Base', 'Var_Base', 'Clinically_Flagged', '0228T_Reads', '5437_NA12878_Reads', '6037_NA12878_Reads']))
        self.assertListEqual(variant_keys, ['Position', 'Ref_Base', 'Var_Base'])
        
    def testMSIParser(self):
        specimens = defaultdict(dict)
        prefixes = []
        variant_keys = []
        control_info=open(path.join(testMSIfile, 'testMSIcontrol'),'rU')
        files = walker(testMSIfile)
        analysis_type='parsers.parse_msi'
        chosen_parser='{}(files, control_info, specimens, prefixes, variant_keys)'.format(analysis_type)
        specimens, prefixes, fieldnames, variant_keys=eval(chosen_parser)  
        self.assertListEqual(sorted(prefixes),sorted(['0228T', '5437_NA12878', '6037_NA12878']))
        self.assertListEqual(sorted(fieldnames), sorted(['0228T', '5437_NA12878', '6037_NA12878', 'Position']))
        self.assertListEqual(variant_keys, ['Position'])
    
      
