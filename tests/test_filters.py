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

files1 = """49_B01_BROv7_NA12878_HA0187_Analysis.txt
49_B01_BROv7_NA12878_HA0187_CNV_Gene_Analysis.txt
49_B01_BROv7_NA12878_HA0187_CNV_Exon_Analysis.txt
49_B01_BROv7_NA12878_HA0187_CNV_QC_Gene_Analysis.txt
49_B01_BROv7_NA12878_HA0187_CNV_QC_Exon_Analysis.txt
49_B01_BROv7_NA12878_HA0187_Pindel_Analysis.txt
49_B01_BROv7_NA12878_HA0187_QC_Analysis.txt
49_B01_BROv7_NA12878_HA0187_SV_Analysis.txt
49_B01_BROv7_NA12878_HA0187_Genotype_Analysis.txt
49_B01_BROv7_NA12878_HA0187_Quality_Analysis.txt
49_B01_BROv7_NA12878_HA0187_MSI_Analysis.txt
49_E04_OPXv4_NA12878_HA0187_Analysis.txt
49_E04_OPXv4_NA12878_HA0187_CNV_Gene_Analysis.txt
49_E04_OPXv4_NA12878_HA0187_CNV_Exon_Analysis.txt
49_E04_OPXv4_NA12878_HA0187_CNV_QC_Gene_Analysis.txt
49_E04_OPXv4_NA12878_HA0187_CNV_QC_Exon_Analysis.txt
49_E04_OPXv4_NA12878_HA0187_Pindel_Analysis.txt
49_E04_OPXv4_NA12878_HA0187_SV_Analysis.txt
49_E04_OPXv4_NA12878_HA0187_Genotype_Analysis.txt
49_E04_OPXv4_NA12878_HA0187_Quality_Analysis.txt
49_E04_OPXv4_NA12878_HA0187_MSI_Analysis.txt
49_C01_OPXv4_HA0187_Analysis.txt
49_C01_OPXv4_HA0187_MSI_Analysis.txt
49_C01_OPXv4_HA0187_CNV_Gene_Analysis.txt
49_C01_OPXv4_HA0187_CNV_Exon_Analysis.txt
49_C01_OPXv4_HA0187_CNV_QC_Gene_Analysis.txt
49_C01_OPXv4_HA0187_CNV_QC_Exon_Analysis.txt
49_C01_OPXv4_HA0187_Pindel_Analysis.txt
49_C01_OPXv4_HA0187_SV_Analysis.txt
49_C01_OPXv4_HA0187_Genotype_Analysis.txt
49_C01_OPXv4_HA0187_Quality_Analysis.txt""".splitlines()


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
        assert keepers == set(['49_B01_BROv7_NA12878_HA0187_Analysis.txt',
                               '49_E04_OPXv4_NA12878_HA0187_Analysis.txt',
                               '49_C01_OPXv4_HA0187_Analysis.txt'])


    def testCNVGeneFileFilter(self):
        keepers = {fn for fn in files1 if f.cnv_gene_analysis(Path('',fn))}
        assert keepers == set(['49_B01_BROv7_NA12878_HA0187_CNV_Gene_Analysis.txt',
                               '49_E04_OPXv4_NA12878_HA0187_CNV_Gene_Analysis.txt',
                               '49_C01_OPXv4_HA0187_CNV_Gene_Analysis.txt'])

    def testCNVExonFileFilter(self):
        keepers = {fn for fn in files1 if f.cnv_exon_analysis(Path('',fn))}
        assert keepers == set(['49_B01_BROv7_NA12878_HA0187_CNV_Exon_Analysis.txt',
                               '49_E04_OPXv4_NA12878_HA0187_CNV_Exon_Analysis.txt',
                               '49_C01_OPXv4_HA0187_CNV_Exon_Analysis.txt'])

    def testPindelFileFilter(self):

        keepers = {fn for fn in files1 if f.pindel_analysis(Path('',fn))}
        assert keepers == set(['49_B01_BROv7_NA12878_HA0187_Pindel_Analysis.txt',
                               '49_E04_OPXv4_NA12878_HA0187_Pindel_Analysis.txt',
                               '49_C01_OPXv4_HA0187_Pindel_Analysis.txt'])

    def testMSIFileFilter(self):

        keepers = {fn for fn in files1 if f.msi_analysis(Path('',fn))}
        assert keepers == set(['49_B01_BROv7_NA12878_HA0187_MSI_Analysis.txt',
                               '49_E04_OPXv4_NA12878_HA0187_MSI_Analysis.txt',
                               '49_C01_OPXv4_HA0187_MSI_Analysis.txt'])

    def testQualityFileFilter(self):

        keepers = {fn for fn in files1 if f.quality_analysis(Path('',fn))}
        assert keepers == set(['49_B01_BROv7_NA12878_HA0187_Quality_Analysis.txt',
                               '49_E04_OPXv4_NA12878_HA0187_Quality_Analysis.txt',
                               '49_C01_OPXv4_HA0187_Quality_Analysis.txt'])
    
    def testGenotypeFileFilter(self):

        keepers = {fn for fn in files1 if f.genotype_analysis(Path('',fn))}
        assert keepers == set(['49_B01_BROv7_NA12878_HA0187_Genotype_Analysis.txt',
                               '49_E04_OPXv4_NA12878_HA0187_Genotype_Analysis.txt',
                               '49_C01_OPXv4_HA0187_Genotype_Analysis.txt'])

