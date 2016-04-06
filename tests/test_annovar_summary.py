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

control ='NA12878-GEN08-HHv1'
NM_dict = {
    'NM_001202435': 'NM_001202435.1',
    'NM_006772': 'NM_006772.1',
    'NM_000038': 'NM_000038.5',
    'NM_007300': 'NM_007300.1',
    'NM_007297': 'NM_007297.2',
    'NM_001015877':'NM_001015877',
    'NM_001123383':'NM_001123383',
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
data4 = {'Gene': 'PHF6(NM_001015877:exon10:c.969-9T>C,NM_032458:exon10:c.969-9T>C),PTEN',
         'Transcripts': '',
         'Variant_Type':'exonic'}
data5 = {'Gene':'BCOR,BCOR(NM_001123383:exon8:c.3503-2A>T,NM_001123384:exon7:c.3449-2A>T,NM_017745:exon8:c.3503-2A>T)',
         'Transcripts':'NM_001123383:exon8:c.3503-2A>T,NM_001123384:exon7:c.3449-2A>T,NM_017745:exon8:c.3503-2A>T)',
         'Variant_Type':'exonic'}
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
        fname = path.join(summary_testfiles, '{}.variant_function').format(control)
        header_ids = {0: 'var_type_1',
                      1: 'gene',
                      7: 'zygosity',
                      12: 'rsid_1',
                      8: 'GATK_Score'}
        variant_idx = [2, 3, 4, 5, 6]
        out = annovar_summary.map_headers(fname, header_ids, variant_idx)
        out = list(out)
        self.assertEqual(len(out), 2546)
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
        data4['Gene'], data4['Transcripts'] = annovar_summary.munge_gene_and_Transcripts(data4, NM_dict)
        data5['Gene'], data5['Transcripts'] = annovar_summary.munge_gene_and_Transcripts(data5, NM_dict)
        #Data 1 gene should be SCN1A
        self.assertTrue(data1['Gene'], 'SCN1A')
        #Data 2 gene should be empyt as the Variant_Type is upstream, which we filter
        self.assertEqual(data2['Gene'], '')
        self.assertEqual(data3['Gene'], 'BRCA1')
        #test that duplicate transcripts are not
        dup_trans='SYNGAP1:NM_006772:exon11:c.1713G>A:p.S571S'
#        self.assertEqual(dup_trans, data2['Transcripts'])

        self.assertEqual(data4['Gene'], 'PHF6,PTEN')
        self.assertEqual(data5['Gene'], 'BCOR')

    def testMungeTranscript(self):
        """
        Return HGVS correct transcript annotations
        Filtered with a preferred transcript list
        NM_006772.1:c.1713G>A
        # """
        data1['c.'], data1['p.'] = annovar_summary.munge_transcript(data1, NM_dict)
        data2['c.'], data2['p.'] = annovar_summary.munge_transcript(data2, NM_dict)
        data3['c.'], data3['p.']  = annovar_summary.munge_transcript(data3, NM_dict)
        data4['c.'], data4['p.']  = annovar_summary.munge_transcript(data4, NM_dict)
        data5['c.'], data5['p.']  = annovar_summary.munge_transcript(data5, NM_dict)
        # #Data 1 gene should be SCN1A
        self.assertEqual(data1['p.'], 'p.A1067T')
        # #Data 2 gene should be empyt as the Variant_Type is upstream, which we filter
        self.assertEqual(data2['c.'], 'NM_006772.1:c.1713G>A')
        #Data 3 p. and c. should have multiple entries
        self.assertIn('p.E1038G',data3['p.']) 
        self.assertIn('p.E991G',data3['p.']) 
        self.assertIn('NM_007297.2:c.2972A>G',data3['c.'])
        self.assertIn('NM_007300.1:c.3113A>G',data3['c.'])

        #4 & 5 should not have a p. but should have a c.
        self.assertIn('NM_001015877:c.969-9T>C',data4['c.'])
        self.assertIn('NM_001123383:c.3503-2A>T',data5['c.'])
        self.assertEqual('',data4['p.'])
        self.assertEqual('',data5['p.'])


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
        self.assertEqual(ref_reads, '-1')
        self.assertEqual(var_reads, '-1')
        self.assertEqual(variant_phred, '')


    def testMungeLJBScores(self):
        """
        Return sift, polyphen and gerp from ljb_all file
        """
        data={'ljb_Scores':'0.012,D,1.0,D,0.81,P,0.092,N,0.999,D,1.355,L,-1.17,T,-1.57,N,0.612,4.883,24.9,0.999,0.919,D,0.022,D,0.529,D,0.706,0,5.51,0.871,0.935,0.826,0.727,16.149'}

        polyphen, sift, mutation_taster, gerp=annovar_summary.munge_ljb_scores(data)
        self.assertEqual(polyphen, '1.0')
        self.assertEqual(sift, '0.012')
        self.assertEqual(mutation_taster,'0.999')
        self.assertEqual(gerp, '5.51')




