"""
Test the masker subcommand script
"""
import os
from os import path
import unittest
import logging
import pprint
import csv
import sys
import json

from munging.subcommands import masker

from __init__ import TestBase
import __init__ as config
log = logging.getLogger(__name__)


analysis_testfiles = path.join(config.datadir, 'analysis_files')
                
class TestMasker(TestBase):
    """
    Test the masker script, which masks by gene name
    """
    def testMaskFileByGene(self):
        """Return only genes listed in the masking dictionary
        """
        data=csv.DictReader(open(path.join(analysis_testfiles,'0228T_CON_OPXv4_INT.SNP_Analysis.txt')), delimiter='\t')
        genes=('BRCA1','BRCA2')
        out_data=masker.mask_file_by_gene(data,genes)
        out_genes=[d['Gene'] for d in out_data]
        self.assertEqual(out_data[0]['Position'],'chr2:12345')
        #2 entries should come out
        self.assertEqual(len(out_data), 2)
        self.assertNotIn('MTHFR', out_genes)
        self.assertIn('BRCA2', out_genes)
