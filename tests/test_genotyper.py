"""
Test the genotyper subcommand script
"""
import os
from os import path
import unittest
import logging
import pprint
import csv
import sys
import json
import pandas as pd
from munging.subcommands import varscan_formatter

from __init__ import TestBase
import __init__ as config
log = logging.getLogger(__name__)

genotyper_testfiles = path.join(config.datadir, 'genotyper')
                
class TestMasker(TestBase):
    """
    Test the varscan_formatter script
    """
    def testFixDeletions(self):
        """Parse deletions/insertions correctly
        """
        cols=['chr','start','stop','ref_base','var_base','comment']
        variants=pd.read_csv(open(path.join(genotyper_testfiles,'Clin_Flagged.txt')), delimiter='\t', header=None, index_col=False, names = cols)
        output=variants.apply(varscan_formatter.fix_deletions, axis = 1)

        self.assertEqual(output.loc[0]['var_base'] ,'INS-9-GGCTCCCCA')
        self.assertEqual(output.loc[1]['var_base'] ,'DEL-1-C')
        self.assertEqual(output.loc[2]['var_base'] ,'DEL-3-GAG')
        self.assertEqual(output.loc[3]['var_base'] ,'G')
        self.assertEqual(output.loc[4]['var_base'] ,'DEL-3-CAT')
        self.assertEqual(output.loc[5]['var_base'] ,'DEL-3-TCT')

        # self.assertNotIn('MTHFR', out_genes)
        # self.assertIn('BRCA2', out_genes)
