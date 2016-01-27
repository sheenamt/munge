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

    def testLaneDetailToSS01(self):
        """Test the lane details are parsed correctly"""
        fcid='HBHW5ADXX'
        ldetail={'FCID':'HBHW5ADXX',
                "Index":'CCAGTTCA', 
                'PlateNumber':'60',
                'Operator':'SF',
                'Well':'D05',
                'SampleProject':'Oncoplex60',
                'Description':'Standard',
                'ControlWell':'E05',
                'Recipe':'OPXv4'}
        
        output=loadlist2samplesheet._lane_detail_to_ss(fcid, ldetail, 1)
        self.assertIn('6036-D05-OPXv4', output)
        #test correct creation of Project name
        self.assertIn('OncoPlex60-OPXv4', output)

    def testGetFlowCellID01(self):
        """Test that only 1 flowcell ID is allowed"""
        reader=csv.DictReader(open(path.join(load_list,'loadlist.csv')))
        #strip whitespace from header names in case tech used wrong template
        reader.fieldnames=[i.strip() for i in reader.fieldnames]
        lane_details = [row for row in reader]
        #SampleSheet.csv needs to be grouped by FCID
        self.assertRaises(ValueError, loadlist2samplesheet._get_flowcell_id, lane_details)

