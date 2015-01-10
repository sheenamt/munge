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

from itertools import ifilter
from collections import namedtuple, defaultdict
from munging import filters 
from munging.utils import Path, walker
from munging.parsers import parse_cnv_exon, parse_cnv_gene, parse_snp, parse_pindel, parse_clin_flagged


from __init__ import TestBase
import __init__ as config
log = logging.getLogger(__name__)

testfiles = path.join(config.datadir, 'analysis_files')

specimens = defaultdict(dict)
annotation = {}
prefixes = []
variant_keys = []
files = ifilter(filters.any_analysis, walker(testfiles))        

class TestParsers(TestBase):
    def setUp(self):
        self.outdir = self.mkoutdir()

    def tearDown(self):
        pass

    def testSNPParser(self):
        """
        Test for correct fieldname parsing 
        """
        analysis_type='parse_snp'
        chosen_parser='{}(files, specimens, annotation, prefixes, variant_keys)'.format(analysis_type)
        specimens, annotation, prefixes, fieldnames, variant_keys=eval(chosen_parser)    
        self.assertListEqual(sorted(prefixes),sorted(['E04_NA12878_Ref|Var','B01_NA12878_Ref|Var','0228T_Ref|Var']))
        self.assertListEqual(fieldnames, ['Position', 'Ref_Base', 'Var_Base', 'Gene', 'Variant_Type',
                                          'Transcripts', 'Clinically_Flagged', 'Cosmic', 'Segdup', 
                                          'Polyphen', 'Sift', 'Mutation_Taster', 'Gerp', 'HiSeq_Freq',
                                          'HiSeq_Count', 'MiSeq_Freq', 'MiSeq_Count', '1000g_ALL', 
                                          'EVS_esp6500_ALL', '1000g_AMR', 'EVS_esp6500_AA', '1000g_EUR',
                                          'EVS_esp6500_EU', '1000g_ASN', '1000g_AFR', 'B01_NA12878_Ref|Var',
                                          'E04_NA12878_Ref|Var', '0228T_Ref|Var'])
        self.assertListEqual(variant_keys, ['Position', 'Ref_Base', 'Var_Base'])
        
#    def testCNVExonParser(self):
 #       analysis_type='parse_cnv_exon'
  #      chosen_parser='{}(files, specimens, annotation, prefixes, variant_keys)'.format(analysis_type)
   #     specimens, annotation, prefixes, fieldnames, variant_keys=eval(chosen_parser)        
    #    #self.assertListEqual(sorted(prefixes),sorted(['E05_Log','0228T_Log']))
     #   self.assertListEqual(variant_keys, ['Position','Gene','Transcripts'])
        #self.assertListEqual(variant_keys, ['Position', 'Ref_Base', 'Var_Base'])
            
    def testCNVGeneParser(self):
        pass
    
    def testQualityParser(self):
        pass
    
    def testPindelParser(self):
        pass
    
      