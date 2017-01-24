"""
Test the genotyper-bam-readcount subcommand script
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
from munging.subcommands import genotyper_bam_readcount

from __init__ import TestBase
import __init__ as config
log = logging.getLogger(__name__)

genotyper_testfiles = path.join(config.datadir, 'genotyper')
                
class TestGenotyper(TestBase):
    """
    Test the genotyper script
    """
    def testFormatIndels(self):
        """Parse deletions/insertions correctly
        BAM requires specially formated variant list where
        deletion == -Variants
        insertion == +Variants
        SNPs are 'start-position', variant
        """
        cols=['chrom','start','stop','Ref_Base','Var_Base','Clinically_Flagged']
        variants=pd.read_csv(open(path.join(genotyper_testfiles,'Clin_Flagged.txt')), delimiter='\t', header=None, index_col=False, names = cols)
        variants['chrom'] = variants['chrom'].astype('str')
        output=variants.apply(genotyper_bam_readcount.format_indels, axis = 1)

        self.assertEqual(output.loc[0]['varscan_variant'] ,'INS-9-GGCTCCCCA')
        self.assertEqual(output.loc[1]['varscan_variant'] ,'DEL-1-C')
        self.assertEqual(output.loc[2]['varscan_variant'] ,'DEL-3-GAG')
        self.assertEqual(output.loc[3]['varscan_variant'] ,'G')
        self.assertEqual(output.loc[4]['varscan_variant'] ,'DEL-3-CAT')
        self.assertEqual(output.loc[5]['varscan_variant'] ,'DEL-3-TCT')

        # self.assertNotIn('MTHFR', out_genes)
        # self.assertIn('BRCA2', out_genes)

    def testParseVarscanLine(self):
        """
        Parse each varscan line into 3 namedtupeles:
        Position
        Reference
        Variant 
        """
        line1='7\t117199640\tT\t16351\t16351\tT:16304:2:26:1:5459:10845:21\tDEL-3-ATC\t\20\t2\t24\t1\t10\t10\tA:3:1:18:1:0:3\tC:10:2:25:1:2:8\tG:5:1:14:1:0:5\tINS-1-A:1:1:26:1:0:1\n'
        
        info1=genotyper.parse_varscan_line(line1)
        self.assertEqual(info1[0][0], '7')
        self.assertEqual(info1[0][1],'117199640')
        self.assertEqual(info1[0][2],'T')
        self.assertEqual(info1[1]['query_variant'][0],'DEL-3-ATC')
