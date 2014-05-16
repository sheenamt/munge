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

from numpy import std, average

from munging.subcommands import summary
from munging.subcommands import annovar_bed_parser
from munging.subcommands import control_parser
from munging.subcommands import qc_variants
from munging.subcommands import quality_metrics
from munging.subcommands import xlsmaker
from munging.subcommands import combined_cnv, combined_pindel, combined_output
from munging.utils import munge_path

from __init__ import TestBase
import __init__ as config
log = logging.getLogger(__name__)


summary_testfiles = path.join(config.datadir, 'summary/annovar_sum')
annovar_testfiles = path.join(config.datadir, 'annovar_bed_parser')
control_testfiles = path.join(config.datadir, 'control_parser')
qc_testfiles = path.join(config.datadir, 'qc_variants')
quality_testfiles = path.join(config.datadir, 'quality_metrics')
varscan_testfiles = path.join(config.datadir, 'varscan')
combined_testfiles = path.join(config.datadir, 'combined')

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
data2 = {'Gene': 'SYNGAP1',
         'Transcripts': 'SYNGAP1:NM_006772:exon11:c.1713G>A:p.S571S',
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
        fname = path.join(summary_testfiles, 'CON-762R.variant_function')
        header_ids = {0: 'var_type_1',
                      1: 'gene',
                      7: 'zygosity',
                      12: 'rsid_1',
                      8: 'GATK_Score'}
        variant_idx = [2, 3, 4, 5, 6]
        out = summary.map_headers(fname, header_ids, variant_idx)
        out = list(out)
        self.assertEqual(len(out), 1373)
        header_keys = set(header_ids.values())
        # confirm that all keys in header_ids are contained in each row of the output
        for pos, data in out:
            self.assertFalse(header_keys - set(data.keys()))

    def testMungeGeneAndTranscripts(self):
        """
        Return modified values of (Gene, Transcripts). Note that
        his depends on 'Variant_Type' provided by munge_variant.
        """
        data1['Gene'], data1['Transcripts'] = summary.munge_gene_and_Transcripts(data1, NM_dict)
        data2['Gene'], data2['Transcripts'] = summary.munge_gene_and_Transcripts(data2, NM_dict)
        data3['Gene'], data3['Transcripts'] = summary.munge_gene_and_Transcripts(data3, NM_dict)

        #Data 1 gene should be SCN1A
        self.assertTrue(data1['Gene'], 'SCN1A')
        #Data 2 gene should be empyt as the Variant_Type is upstream, which we filter
        self.assertEqual(data2['Gene'], '')

    def testMungeTranscript(self):
        """
        Return HGVS correct transcript annotations
        Filtered with a preferred transcript list
        NM_006772.1:c.1713G>A
        # """
        data1['c.'], data1['p.'] = summary.munge_transcript(data1, NM_dict)
        data2['c.'], data2['p.'] = summary.munge_transcript(data2, NM_dict)
        data3['c.'], data3['p.'] = summary.munge_transcript(data3, NM_dict)

        # #Data 1 gene should be SCN1A
        self.assertEqual(data1['p.'], 'p.A1067T')
        # #Data 2 gene should be empyt as the Variant_Type is upstream, which we filter
        self.assertEqual(data2['c.'], 'NM_006772.1:c.1713G>A')
        #Data 3 p. and c. should have multiple entries
        self.assertEqual(data3['p.'], 'p.E1038G p.E991G')
        self.assertEqual(data3['c.'], 'NM_007300.1:c.3113A>G NM_007297.2:c.2972A>G')

    def testGetAlleleFreq(self):
        """
        Return allele frequency of var_reads/ref_reads
        """
        freq1 = summary.get_allele_freq(data1)
        freq2 = summary.get_allele_freq(data2)
        self.assertEqual(freq1, 'NA')
        self.assertEqual(freq2, '0.10')

    def testGetReads(self):
        """
        Return Ref_Reads, Var_Reads, Variant_Phred
        """
        data='1/1:255:289:289:18:267:92.39%:3.3502E-142:24:30:7:11:48:219'
        ref_reads, var_reads, variant_phred = summary.get_reads(data)
        self.assertEqual(ref_reads, '18')
        self.assertEqual(var_reads, '267')
        self.assertEqual(variant_phred, '30')



class TestAnnovarBedParser(TestBase):
    """
    Test the annovar_bed_parser which is used to filter the vcf/varscan/varscanSNP
    files to include only the regions specified in the bed file for annotation
    """

    def testAnnovarBedParserCoords(self):
        """
        Takes in row of file, return chr, start, stop
        """
        bedfname = open(path.join(annovar_testfiles, 'test.bed'))
        bedinfo = list(csv.reader(bedfname, delimiter='\t'))
        for row in bedinfo:
            bedoutput = annovar_bed_parser.coords(row)
            self.assertEqual(len(bedoutput), 3)


class TestControlParser(TestBase):
    """
    Test the control_parser which is used to check the control sample against
    SNVs found in NA12878 from Complete Genomics, 1000G and our samples
    """

    def testControlParserMatch(self):
        """
        Make a list and keep count of variants found in both the qc file
        and the run output for LMG/OPX-240 sample
        Matches if chr, start are the same
        control[chr] = run[chr] and control[start] = run[start]
        """
        controlfname = open(path.join(control_testfiles, 'OncoPlex_qc_variants_v3.txt'))
        controlinfo = list(csv.reader(controlfname, delimiter='\t'))
        runfname = open(path.join(control_testfiles, 'OPX-240_Analysis.txt'))
        runinfo = list(csv.reader(runfname, delimiter='\t'))
        output, count = control_parser.match(controlinfo, runinfo)
        #Count and output length should be qual
        self.assertEqual(len(output), count)
        #The second entry of the second line should be MTHFR:NM_005957:exon8:c.1286A>C:p.E429A,
        self.assertEqual(output[1][1], 'MTHFR:NM_005957:exon8:c.1286A>C:p.E429A,')


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
        pipefname = open(path.join(qc_testfiles, 'OPX-240.exonic_variant_function'))
        pipe = list(csv.reader(pipefname, delimiter="\t"))
        kgfname = open(path.join(qc_testfiles, 'NA12878.1000g.hg19.exonic_variant_function'))
        kg = list(csv.reader(kgfname, delimiter="\t"))
        cgfname = open(path.join(qc_testfiles, 'NA12878.CG.hg19.exonic_variant_function'))
        cg = list(csv.reader(cgfname, delimiter="\t"))
        output = qc_variants.match(pipe, kg, cg)
        #There should be 3 lines that match
        self.assertEqual(len(output), 3)
        #The second entry on the second line should be SYNGAP1:NM_006772:exon11:c.1713G>A:p.S571S
        self.assertEqual(str(output[1][1]), "SYNGAP1:NM_006772:exon11:c.1713G>A:p.S571S,")


class TestQualityMetrics(TestBase):
    """
    Test the quality_metrics script with reports the average read depth,
    standard devation given the CNV_bins file and picard metrics file
    """

    def testGetValues(self):
        """
        Get the standard deviation from the tumor.rd.ori column in the CNV file
        """
        cnvfname = open(path.join(quality_testfiles, 'OPX-240_CNV_bins.txt'))
        cnv_info = csv.reader(cnvfname, delimiter="\t")
        cnv_info.next()
        a = quality_metrics.get_values(cnv_info)
        stdev = ["Standard deviation of read depth: %0.2f " % (a).std()]
        ave = ["Average read depth: %0.2f " % average(a)]
        self.assertEqual(stdev, ['Standard deviation of read depth: 120.25 '])
        self.assertEqual(ave, ['Average read depth: 312.25 '])


class TestXlsmaker(TestBase):
    """
    Test the xlsmaker which combines all the analysis files
    into one workbook and renames sheets for clarity in sign out
    """

    def testfloatifpossible(self):
        """
        Convert integers to float instead of string where applicable.
        """
        test01 = 'Gene'
        test02 = '12'
        self.assertTrue(xlsmaker.float_if_possible(test01), 'Gene')
        self.assertTrue(xlsmaker.float_if_possible(test02), '12.0')

    def testProcessFiles(self):
        """
        Rename the analysis files for workbook
        """
        tab = '10_SNP_Indel'
        filetype = 'Analysis'
        files = []
        files.append(path.join(summary_testfiles, 'LMG-098A_Analysis.txt'))
        files.append(path.join(summary_testfiles, 'LMG-098A_Quality_Analysis.txt'))
        data, fname = xlsmaker.process_files(files, tab, filetype)
        self.assertEqual(data, '10_SNP_Indel')
        self.assertEqual(fname, 'testfiles/summary/annovar_sum/LMG-098A_Analysis.txt')
        self.assertNotEqual(fname, 'testfiles/summary/annovar_sum/LMG-098A_Quality_Analysis.txt')



