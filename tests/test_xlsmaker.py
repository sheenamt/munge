"""
Test the subcommand scripts
"""
import os
from os import path
import unittest
import logging
import pprint
import csv
import sys
import json

from munging.subcommands import xlsmaker

from __init__ import TestBase
import __init__ as config
log = logging.getLogger(__name__)

summary_testfiles = path.join(config.datadir, 'annovar_summary')
annovar_testfiles = path.join(config.datadir, 'annovar_bed_parser')
control_testfiles = path.join(config.datadir, 'control_parser')
qc_testfiles = path.join(config.datadir, 'qc_variants')
quality_testfiles = path.join(config.datadir, 'quality_metrics')
analysis_testfiles = path.join(config.datadir, 'analysis_files')
load_list = path.join(config.datadir, 'loadlist')

control ='5437_E05_OPXv4_NA12878_MA0013'
NM_dict = {
    'NM_001202435': 'NM_001202435.1',
    'NM_006772': 'NM_006772.1',
    'NM_000038': 'NM_000038.5',
    'NM_007300': 'NM_007300.1',
    'NM_007297': 'NM_007297.2'
}
data1 = {'Gene': 'SCN1A',
         'Transcripts': 'SCN1A:NM_001202435:exon18:c.3199G>A:p.A1067T,',
         'Variant_Type': '',
         'Var_Reads': '-1', 'Ref_Reads': '-1'}
#Duplicate transcript entry to test parsing of duplicates
data2 = {'Gene': 'SYNGAP1',
         'Transcripts': 'SYNGAP1:NM_006772:exon11:c.1713G>A:p.S571S,SYNGAP1:NM_006772:exon11:c.1713G>A:p.S571S,SYNGAP1:NM_006772:exon11:c.1713G>A:p.S571S',
         'Variant_Type': 'upstream',
         'Var_Reads': '10', 'Ref_Reads': '90'}
data3 = {'Gene': 'BRCA1',
         'Transcripts': 'BRCA1:NM_007300:exon10:c.3113A>G:p.E1038G,BRCA1:NM_007297:exon9:c.2972A>G:p.E991G,BRCA1:NM_007294:exon10:c.3113A>G:p.E1038G',
         'Variant_Type': ''}


class TestXlsmaker(TestBase):
    """
    Test the xlsmaker which combines all the analysis files
    into one workbook and renames sheets for clarity in sign out
    """

    def testfloatifpossible(self):
        """
        Convert integers to float instead of string where applicable.
        """
        test01 = 'Gene'
        test02 = '12'
        self.assertTrue(xlsmaker.float_if_possible(test01), 'Gene')
        self.assertTrue(xlsmaker.float_if_possible(test02), '12.0')

    def testProcessFiles(self):
        """
        Rename the analysis files for workbook
        """
        tab = '10_SNP_Indel'
        filetype = 'Analysis'
        files = []
        files.append(path.join(summary_testfiles, '{}.SNP_Analysis.txt'.format(control)))
        files.append(path.join(summary_testfiles, '{}.Quality_Analysis.txt'.format(control)))
        data, fname = xlsmaker.process_files(files, tab, filetype)
        self.assertEqual(data, '10_SNP_Indel')
        self.assertEqual(fname, 'testfiles/annovar_summary/{}.SNP_Analysis.txt'.format(control))
        self.assertNotEqual(fname, 'testfiles/annovar_summary/{}.Quality_Analysis.txt'.format(control))

