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

from collections import namedtuple
from munging import filters as f
from munging.utils import Path

from __init__ import TestBase
import __init__ as config
log = logging.getLogger(__name__)

files1 = """LMG-240-1-10_v2_1_Analysis.txt
LMG-240-1-10_v2_1_CNV_Gene_Analysis.txt
LMG-240-1-10_v2_1_CNV_Exon_Analysis.txt
LMG-240-1-10_v2_1_CNV_QC_Analysis.txt
LMG-240-1-10_v2_1_Pindel_Analysis.txt
LMG-240-1-10_v2_1_QC_Analysis.txt
LMG-240-1-10_v2_1_SV_Analysis.txt
LMG-240-1-10_v2_1_Genotype_Analysis.txt
LMG-240-1-10_v2_1_Quality_Analysis.txt
OPX-0124T_Analysis.txt
OPX-0124T_CNV_Gene_Analysis.txt
OPX-0124T_CNV_Exon_Analysis.txt
OPX-0124T_CNV_QC_Analysis.txt
OPX-0124T_Pindel_Analysis.txt
OPX-0124T_SV_Analysis.txt
OPX-0124T_Genotype_Analysis.txt
OPX-0124T_Quality_Analysis.txt
OPX-0129T_Analysis.txt
OPX-0129T_CNV_Gene_Analysis.txt
OPX-0129T_CNV_Exon_Analysis.txt
OPX-0129T_CNV_QC_Analysis.txt
OPX-0129T_Pindel_Analysis.txt
OPX-0129T_SV_Analysis.txt
OPX-0129T_Genotype_Analysis.txt
OPX-0129T_Quality_Analysis.txt""".splitlines()


class TestFilters(TestBase):


    def setUp(self):
        self.outdir = self.mkoutdir()

    def tearDown(self):
        pass

    def testAnyAnalysisFilter(self):
        assert f.any_analysis(Path('','Analysis.txt')) is True
        assert f.any_analysis(Path('','Analysis.csv')) is True
        assert f.any_analysis(Path('','Analsis.csv')) is False
        assert f.any_analysis(Path('','Analysis.xls')) is False
        assert f.any_analysis(Path('','Analysis.xlsx')) is False


    def testOnlyAnalysisFilter(self):

        keepers = {fn for fn in files1 if f.only_analysis(Path('',fn))}
        assert keepers == set(['LMG-240-1-10_v2_1_Analysis.txt',
                               'OPX-0124T_Analysis.txt',
                               'OPX-0129T_Analysis.txt'])


    def testCNVGeneFileFilter(self):
        keepers = {fn for fn in files1 if f.cnv_gene_analysis(Path('',fn))}
        assert keepers == set(['LMG-240-1-10_v2_1_CNV_Gene_Analysis.txt',
                               'OPX-0124T_CNV_Gene_Analysis.txt',
                               'OPX-0129T_CNV_Gene_Analysis.txt'])

    def testCNVExonFileFilter(self):
        keepers = {fn for fn in files1 if f.cnv_exon_analysis(Path('',fn))}
        assert keepers == set(['LMG-240-1-10_v2_1_CNV_Exon_Analysis.txt',
                               'OPX-0124T_CNV_Exon_Analysis.txt',
                               'OPX-0129T_CNV_Exon_Analysis.txt'])

    def testPindelFileFilter(self):

        keepers = {fn for fn in files1 if f.pindel_analysis(Path('',fn))}
        assert keepers == set(['LMG-240-1-10_v2_1_Pindel_Analysis.txt',
                               'OPX-0124T_Pindel_Analysis.txt',
                               'OPX-0129T_Pindel_Analysis.txt'])

