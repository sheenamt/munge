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

from munging.subcommands import annovar_bed_parser

from __init__ import TestBase
import __init__ as config
log = logging.getLogger(__name__)

annovar_testfiles = path.join(config.datadir, 'annovar_bed_parser')


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

