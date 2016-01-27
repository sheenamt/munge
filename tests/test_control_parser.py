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

from munging.subcommands import control_parser

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

class TestControlParser(TestBase):
    """
    Test the control_parser which is used to check the control sample against
    SNVs found in NA12878 from Complete Genomics, 1000G and our samples
    """

    def testControlParserMatch(self):
        """
        Make a list and keep count of variants found in both the qc file
        and the run output for LMG/OPX-240 sample
        Matches if chr, start are the same
        control[chr] = run[chr] and control[start] = run[start]
        """
        controlfname = open(path.join(control_testfiles, 'OncoPlex_qc_variants_v4.txt'))
        controlinfo = list(csv.reader(controlfname, delimiter='\t'))
        runfname = open(path.join(analysis_testfiles,'{}.SNP_Analysis.txt').format(control))
        runinfo = list(csv.reader(runfname, delimiter='\t'))
        output, count = control_parser.match(controlinfo, runinfo)
        #Count and output length should be qual
        self.assertEqual(len(output), count)
        #The second entry of the second line should be MTHFR:NM_005957:exon8:c.1286A>C:p.E429A,
        self.assertEqual(output[0][1], 'MTHFR:NM_005957:exon8:c.1305C>T:p.F435F,')


