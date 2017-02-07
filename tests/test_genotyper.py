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
from munging.subcommands import genotyper

from __init__ import TestBase
import __init__ as config
log = logging.getLogger(__name__)

genotyper_testfiles = path.join(config.datadir, 'genotyper')
                
class TestMasker(TestBase):
    """
    Test the varscan_formatter script
    """
    def testFormatIndels(self):
        """Parse deletions/insertions correctly
         Deletions are 'start-position'-1, DEL-len(reference)-reference
         Insertions are 'start-position', INS-len(variant)-variant
         SNPs are 'start-position', variant
        """
        cols=['chrom','start','stop','Ref_Base','Var_Base','Clinically_Flagged']
        variants=pd.read_csv(open(path.join(genotyper_testfiles,'Clin_Flagged.txt')), delimiter='\t', header=None, index_col=False, names = cols)
        variants['chrom'] = variants['chrom'].astype('str')
        output=variants.apply(genotyper.format_indels, axis = 1)

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
        Variant = namedtuple('Variant', ['base','reads','strands','avg_qual','map_qual','plus_reads','minus_reads'])
        Position = namedtuple('Position', ['chrom','position','ref_base', 'depth'])
        Reference = namedtuple('Reference', ['base','reads','strands','avg_qual','map_qual','plus_reads','minus_reads', 'misc'
        """
        line1='7\t117199640\tT\t16351\t16301\tT:16280:2:26:1:5459:10845:21\tDEL-3-ATC\t\20\t2\t24\t1\t10\t10\tA:3:1:18:1:0:3\tC:10:2:25:1:2:8\tG:5:1:14:1:0:5\tINS-1-A:1:1:26:1:0:1\n'
        
        info1=genotyper.parse_varscan_line(line1)
        variant1='DEL-3-ATC'
        variant2='INS-1-A'
        variant3='G'
        variant4='A'
        variant5='C'
        #Chrom
        self.assertEqual(info1[0][0], '7')
        #position
        self.assertEqual(info1[0][1],'117199640')
        #refbase
        self.assertEqual(info1[0][2],'T')
        #depth with qfilter (set in readcounts call)
        self.assertEqual(info1[0][3],'16301')
        #check that multiple variants from one position are found
        for variant in info1[1]['variants']:
            if variant[0]==variant1:
                self.assertEqual(variant[0],'DEL-3-ATC')
                self.assertEqual(variant[1],'20')
            if variant[0]==variant2:
                self.assertEqual(variant[0],'INS-1-A')
                self.assertEqual(variant[1],'1')
            if variant[0]==variant3:
                self.assertEqual(variant[0],'G')
                self.assertEqual(variant[1],'5')
            if variant[0]==variant4:
                self.assertEqual(variant[0],'A')
                self.assertEqual(variant[1],'3')
            if variant[0]==variant5:
                self.assertEqual(variant[0],'C')
                self.assertEqual(variant[1],'10')

