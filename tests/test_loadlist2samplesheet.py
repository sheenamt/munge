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

from munging.subcommands import loadlist2samplesheet
from munging.subcommands.masker import MASK_CODES

from __init__ import TestBase
import __init__ as config
log = logging.getLogger(__name__)

load_list = path.join(config.datadir, 'loadlist')

class TestLoadListtoSampleSheet(TestBase):
    """
    Test the script with produces the demux sample sheet
    """
    genes = ['GATA2', 'BRCA1', 'BRCA2']
    fcid='HBHW5ADXX'
    ldetail={'FCID':'HBHW5ADXX',
             "Index":'CCAGTTCA', 
             'PlateNumber':'60',
             'Operator':'SF',
             'Well':'D05',
             'SampleProject':'EpiPlex60',
             'Description':'Standard',
             'ControlWell':'E05',
             'Recipe':'EPIv2',
             'CPT':'EPIPX1',
             'Genes':['GATA2',]}
        

    def testLaneDetailToSS01(self):
        """Test the lane details are parsed correctly"""
        output=loadlist2samplesheet._lane_detail_to_ss(self.fcid, self.ldetail, 1, self.genes)
        self.assertIn('6036-D05-EPIv2', output)
        #test correct creation of Project name
        self.assertIn('EpiPlex60-EPIv2', output)

    def testGetFlowCellID01(self):
        """Test that only 1 flowcell ID is allowed"""
        reader=csv.DictReader(open(path.join(load_list,'loadlist.csv')))
        #strip whitespace from header names in case tech used wrong template
        reader.fieldnames=[i.strip() for i in reader.fieldnames]
        lane_details = [row for row in reader]
        #SampleSheet.csv needs to be grouped by FCID
        self.assertRaises(ValueError, loadlist2samplesheet._get_flowcell_id, lane_details)


    def testCPTandGeneToList(self):
        """Test the parsing of the cpt and gene codes"""
        #Valid mask code, should return mask code
        self.assertEqual(loadlist2samplesheet.cpt_and_gene_to_list(self.ldetail, self.genes),'EPIPX1')
        #CPT not a valid masking, so gene list should be empty
        ldetail2={'CPT':'BROCA',
                  'Recipe': 'BROv10',
                  'Genes':''}
        self.assertEqual(loadlist2samplesheet.cpt_and_gene_to_list(ldetail2, self.genes), '')
        #If empty CPT, provide gene list
        ldetail3={'CPT':'',
                  'Recipe':'EPIv2',
                  'Genes':['GATA2','BRCA1']}
        self.assertEqual(loadlist2samplesheet.cpt_and_gene_to_list(ldetail3, self.genes), ['GATA2','BRCA1'])
