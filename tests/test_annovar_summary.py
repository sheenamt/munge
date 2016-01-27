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

from munging.subcommands import annovar_summary

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

class TestSummary(TestBase):
    """
    Test the summary script with creates the Analysis.txt file from
    annotations created with annovar.
    """

    def setUp(self):
        self.outdir = self.mkoutdir()

    def testMapHeaders(self):
        """
        Gets header(s) and info from each file
        """
        fname = path.join(summary_testfiles, '{}.gatk.variant_function').format(control)
        header_ids = {0: 'var_type_1',
                      1: 'gene',
                      7: 'zygosity',
                      12: 'rsid_1',
                      8: 'GATK_Score'}
        variant_idx = [2, 3, 4, 5, 6]
        out = annovar_summary.map_headers(fname, header_ids, variant_idx)
        out = list(out)
        self.assertEqual(len(out), 1265)
        header_keys = set(header_ids.values())
        # confirm that all keys in header_ids are contained in each row of the output
        for pos, data in out:
            self.assertFalse(header_keys - set(data.keys()))

    def testMungeGeneAndTranscripts(self):
        """
        Return modified values of (Gene, Transcripts). Note that
        his depends on 'Variant_Type' provided by munge_variant.
        """
        data1['Gene'], data1['Transcripts'] = annovar_summary.munge_gene_and_Transcripts(data1, NM_dict)
        data2['Gene'], data2['Transcripts'] = annovar_summary.munge_gene_and_Transcripts(data2, NM_dict)
        data3['Gene'], data3['Transcripts'] = annovar_summary.munge_gene_and_Transcripts(data3, NM_dict)

        #Data 1 gene should be SCN1A
        self.assertTrue(data1['Gene'], 'SCN1A')
        #Data 2 gene should be empyt as the Variant_Type is upstream, which we filter
        self.assertEqual(data2['Gene'], '')
        self.assertEqual(data3['Gene'], 'BRCA1')
        #test that duplicate transcripts are not
        dup_trans='SYNGAP1:NM_006772:exon11:c.1713G>A:p.S571S'
        self.assertEqual(dup_trans, data2['Transcripts'])


    def testMungeTranscript(self):
        """
        Return HGVS correct transcript annotations
        Filtered with a preferred transcript list
        NM_006772.1:c.1713G>A
        # """
        data1['c.'], data1['p.'] = annovar_summary.munge_transcript(data1, NM_dict)
        data2['c.'], data2['p.'] = annovar_summary.munge_transcript(data2, NM_dict)
        data3['c.'], data3['p.']  = annovar_summary.munge_transcript(data3, NM_dict)

        # #Data 1 gene should be SCN1A
        self.assertEqual(data1['p.'], 'p.A1067T')
        # #Data 2 gene should be empyt as the Variant_Type is upstream, which we filter
        self.assertEqual(data2['c.'], 'NM_006772.1:c.1713G>A')
        #Data 3 p. and c. should have multiple entries
        self.assertEqual('p.E991G p.E1038G',data3['p.']) #, 'p.E1038G p.E991G')
        self.assertEqual(data3['c.'], 'NM_007297.2:c.2972A>G NM_007300.1:c.3113A>G')

    def testGetAlleleFreq(self):
        """
        Return allele frequency of var_reads/ref_reads
        """
        freq1 = annovar_summary.get_allele_freq(data1)
        freq2 = annovar_summary.get_allele_freq(data2)
        self.assertEqual(freq1, 'NA')
        self.assertEqual(freq2, '0.10')

    def testGetReads(self):
        """
        Return Ref_Reads, Var_Reads, Variant_Phred
        """
        varscan_header='GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR'
        data1='1/1:255:289:289:18:267:92.39%:3.3502E-142:24:30:7:11:48:219'
        ref_reads, var_reads, variant_phred = annovar_summary.get_reads(varscan_header,data1)
        self.assertEqual(ref_reads, '18')
        self.assertEqual(var_reads, '267')
        self.assertEqual(variant_phred, '30')

        data2='1/1:0,2:2:6:65,6,0'
        gatk_header='GT:AD:DP:GQ:PL'
        ref_reads, var_reads, variant_phred = annovar_summary.get_reads(gatk_header,data2)
        self.assertEqual(ref_reads, '0')
        self.assertEqual(var_reads, '2')
        self.assertEqual(variant_phred, '')


    def testMungeLJBScores(self):
        """
        Return sift, polyphen and gerp from ljb_all file
        """
        data={'ljb_Scores':'0.18,T,0.462,P,0.065,B,0.000,D,0.000,P,1.375,L,-0.78,T,-1.076,T,0.000,T,0.397,2.791,15.29,4.74,1.898,3.792,13.742'}

        polyphen, sift, gerp, mutation_taster=annovar_summary.munge_ljb_scores(data)
        self.assertEqual(polyphen, '0.462')
        self.assertEqual(sift, '0.18')
        self.assertEqual(gerp, '4.74')
        self.assertEqual(mutation_taster,'0.000')



