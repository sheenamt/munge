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

        self.assertEqual(output.loc[0]['readcount_variant'] ,'+GGCTCCCCA')
        self.assertEqual(output.loc[1]['readcount_variant'] ,'-C')
        self.assertEqual(output.loc[2]['readcount_variant'] ,'-GAG')
        self.assertEqual(output.loc[3]['readcount_variant'] ,'G')
        self.assertEqual(output.loc[4]['readcount_variant'] ,'-CAT')
        self.assertEqual(output.loc[5]['readcount_variant'] ,'-TCT')

    def testParseReadCountLine(self):
        """
        Parse each varscan line into 3 namedtupeles:
        Position
        Reference
        Variant 
        """
        line1='2\t234668879\tC\t278\t=:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00\tA:2:44.50:10.50:33.00:1:1:0.43:0.07:48.00:1:0.05:94.00:0.18\tC:272:56.75:27.72:36.11:117:155:0.54:0.02:12.65:117:0.44:100.06:0.49\tG:2:70.00:8.50:37.00:1:1:0.67:0.04:20.00:1:0.30:101.00:0.33\tT:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00\tN:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00\t+AT:80:58.97:0.00:36.28:33:47:0.62:0.03:9.49:33:0.47:99.97:0.48\t+ATAT:7:60.00:0.00:37.00:3:4:0.60:0.04:2.43:3:0.30:101.00:0.31\t-CA:2:29.00:0.00:29.00:1:1:0.69:0.05:24.00:1:0.50:84.50:0.30\n'
        
        info=genotyper_bam_readcount.parse_readcount_line(line1)
        chrom = info[0][0]
        pos_start = int(info[0][1])
        depth = info[0][3]
        Reference_Reads=info[2]['reference'][0][1]
        variant1='+ATAT'
        variant2='-CA'
        self.assertEqual(chrom, '2')
        self.assertEqual(pos_start,234668879)
        self.assertEqual(depth,'278')
        self.assertEqual(Reference_Reads,'272')
        for variant in info[1]['variants']:
            if variant[0]==variant1:
                self.assertEqual(variant[0],'+ATAT')
                self.assertEqual(variant[1],'7')
            if variant[0]==variant2:
                self.assertEqual(variant[0],'-CA')
                self.assertEqual(variant[1],'2')
