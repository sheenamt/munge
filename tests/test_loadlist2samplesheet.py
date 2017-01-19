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

from __init__ import TestBase
import __init__ as config
log = logging.getLogger(__name__)

load_list = path.join(config.datadir, 'loadlist')

class TestLoadListtoSampleSheet(TestBase):
    """
    Test the script with produces the demux sample sheet
    """
    fcid='HBHW5ADXX'
    ldetail1={'FCID':'HBHW5ADXX',
             "Index":'CCAGTTCA', 
             'PlateNumber':'60',
             'Operator':'SF',
             'Well':'D05',
             'SampleProject':'EpiPlex60',
             'Description':'Standard',
             'ControlWell':'E05',
             'Recipe':'EPIv2'}

    def testLaneDetailToSS01(self):
        """Test the lane details are parsed correctly"""
        output=loadlist2samplesheet._lane_detail_to_ss(self.fcid, self.ldetail1, 1)
        self.assertIn('6036-D05-EPIv2', output)
        #test correct creation of Project name
        self.assertIn('EpiPlex60-EPIv2', output)

    def testGetFlowCellID01(self):
        """Test that only 1 flowcell ID is allowed"""
        reader=csv.DictReader(open(path.join(load_list,'bad-loadlist1.csv')))
        #strip whitespace from header names in case tech used wrong template
        reader.fieldnames=[i.strip() for i in reader.fieldnames]
        lane_details = [row for row in reader]
        #SampleSheet.csv needs to be grouped by FCID
        self.assertRaises(ValueError, loadlist2samplesheet._get_flowcell_id, lane_details)

    def testControlCheck(self):
        """Test that 1 control required for each assay"""
        reader=csv.DictReader(open(path.join(load_list,'bad-loadlist2.csv')))
        #strip whitespace from header names in case tech used wrong template
        reader.fieldnames=[i.strip() for i in reader.fieldnames]
        lane_details = [row for row in reader]
        #SampleSheet.csv needs 1 NA12878 per assay
        self.assertRaises(ValueError, loadlist2samplesheet.check_control_vs_number_assays, lane_details)

    def testSortByWell(self):
        """Test that the input is sorted on Well, which drives the name created later"""
        reader=csv.DictReader(open(path.join(load_list,'loadlist.csv')))
        #strip whitespace from header names in case tech used wrong template
        reader.fieldnames=[i.strip() for i in reader.fieldnames]
        lane_details = [row for row in reader]
        sorted_lane_details = sorted(lane_details, key=lambda d: d['Well'])
        #SampleSheet.csv needs to be grouped by FCID
        input_wells = ['A01', 'D01', 'B01', 'C01', 'E01']
        not_sorted_wells = [f['Well'] for f in lane_details]
        output_wells = [f['Well'] for f in sorted_lane_details]
        expected_wells = ['A01', 'B01', 'C01', 'D01', 'E01']
        self.assertListEqual(input_wells, not_sorted_wells)
        self.assertListEqual(output_wells, expected_wells)
        
