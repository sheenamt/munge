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

from munging.subcommands import qc_variants

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


class TestQCVariants(TestBase):
    """
    Test the qc_variants script which finds the intersection of
    (1000G, Complete Genomes, LMG/OPX-240 output)
    to create versioned assay specific qc file
    """

    def testQCVariantsMatch(self):
        """
        Return list of variants found in all three input files
        (1000G, Complete Genomes, LMG/OPX-240 output)
        Matches if chrm, start, stop, ref_base, and var_base are the same.
        """
        pipefname = open(path.join(qc_testfiles, '{}.merged.exonic_variant_function').format(control))
        pipe = list(csv.reader(pipefname, delimiter="\t"))
        kgfname = open(path.join(qc_testfiles, 'NA12878.1000g.hg19.exonic_variant_function'))
        kg = list(csv.reader(kgfname, delimiter="\t"))
        cgfname = open(path.join(qc_testfiles, 'NA12878.CG.hg19.exonic_variant_function'))
        cg = list(csv.reader(cgfname, delimiter="\t"))
        output = qc_variants.match(pipe, kg, cg)
        #There should be 3 lines that match
        self.assertEqual(len(output), 3)
        #The second entry on the second line should be SCN8A:NM_001177984:exon5:c.576C>T:p.D192D,SCN8A:NM_014191:exon5:c.576C>T:p.D192D
        self.assertEqual(str(output[1][1]), "SCN8A:NM_001177984:exon5:c.576C>T:p.D192D,SCN8A:NM_014191:exon5:c.576C>T:p.D192D,")

