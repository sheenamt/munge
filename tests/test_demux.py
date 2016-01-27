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

from munging.subcommands import demux

from __init__ import TestBase
import __init__ as config
log = logging.getLogger(__name__)


class TestDemux(TestBase):
    """
    Test the flowcell directory parser
    """

    def testParseFlowcellDir(self):
        """Test the flowcell info is parsed correctly"""
        fcid1='/media/NGS-Data/NextSeq/150203_NS500359_0004_AH03JVAFXX'
        fcid2='/home/illumina/miseq/141023_M00829_0003_000000000-AAC7L'
        fcid3='/home/illumina/hiseq/140902_D00180_0185_AHAC3AADXX'
        demux.parse_flowcell_dir(fcid1)
        demux.parse_flowcell_dir(fcid2)
        demux.parse_flowcell_dir(fcid3)

        self.assertDictEqual(dict(demux.parse_flowcell_dir(fcid1)), 
                             {'flowcell_dir':'150203_NS500359_0004_AH03JVAFXX',
                              'run_date': '150203',
                              'Illumina_id': 'NS500359',
                              'machine_run': '0004',
                              'machine_side': 'A',
                              'flowcell_id': 'H03JVAFXX',
                              'drive': '/media/NGS-Data/NextSeq'})
        self.assertDictEqual(dict(demux.parse_flowcell_dir(fcid2)),
                             {'flowcell_dir': '141023_M00829_0003_000000000-AAC7L',
                              'run_date': '141023',
                              'Illumina_id': 'M00829',
                              'machine_run': '0003',
                              'machine_side': '-',
                              'flowcell_id': 'AAC7L',
                              'drive': '/home/illumina/miseq'})
        self.assertDictEqual(dict(demux.parse_flowcell_dir(fcid3)),
                             {'flowcell_dir': '140902_D00180_0185_AHAC3AADXX',
                              'run_date': '140902',
                              'Illumina_id': 'D00180',
                              'machine_run': '0185',
                              'machine_side': 'A',
                              'flowcell_id': 'HAC3AADXX',
                              'drive': '/home/illumina/hiseq'})
