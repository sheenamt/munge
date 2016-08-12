"""
Test the hgvs_check subcommand script
"""

import subprocess
import filecmp
import logging
import os

from munging.subcommands import hgvs_check

from __init__ import TestBase
import __init__ as config
log = logging.getLogger(__name__)


analysis_testfiles = os.path.join(config.datadir, 'analysis_files')
hgvs_testfiles = os.path.join(config.datadir, 'hgvs_files')
control_sample = '6037_E05_OPXv4_NA12878_HA0201'


class TestHGVSCheck(TestBase):
    """
    Test the hgvs_check script, which uploads HGVS annotations to mutalyzer.nl's
    batch web service.
    """
    def setUp(self):
        url = hgvs_check.MUTALYZER_URL
        self.conn = hgvs_check.get_mutalyzer_connection(url)

    def testGetMutalyzerConnection(self):
        """Should be able to connect to mutalyzer.nl
        """
        r = self.conn.ping()
        self.assertEqual(r, "pong")

    def testQueryMutalyzer(self):
        hgvs_test_vars = ["NM_017668.2:c.*883_*884insA",
                          "NM_006019.3:c.2274_2275insGGCCTG",
                          "NM_006622.3:c.*551_*552insT",
                          "NM_002451.3:c.*250delA"]
        r = hgvs_check.query_mutalyzer(self.conn, hgvs_test_vars)
        r = r[1:-1]  # exclude header and additional last row
        fixed_hgvs = [l.split("\t")[6] for l in r]
        known_output = ['c.*898dup', 'c.2285_2290dup', 'c.*551dup', 'c.*258del']
        self.assertEqual(fixed_hgvs, known_output)

    def testCleanProteinField(self):
        hgvs_test_protein_vars = ['p.(=)', 'p.(Leu762_Gly763dup)', 'p.Pro197_Pro198dup']
        cleaned_output = [hgvs_check.clean_protein_field(x) for x in hgvs_test_protein_vars]
        known_output = ['', 'p.L762_G763dup', 'p.P197_P198dup']
        self.assertEqual(cleaned_output, known_output)

    def testHGVSCheck(self):
        """
        Command used to generate expected output:
        ./munge hgvs_check \
            --update --output-error-column \
            -i testfiles/analysis_files/6037_E05_OPXv4_NA12878_HA0201.SNP_Analysis.txt \
            -o testfiles/hgvs_files/6037_E05_OPXv4_NA12878_HA0201.SNP_Analysis.hgvs_fixed.txt
        """
        infile = os.path.join(analysis_testfiles, '{}.SNP_Analysis.txt').format(control_sample)
        outfile = os.path.join(self.mkoutdir(), 'temp.txt')

        expected = os.path.join(hgvs_testfiles,
                                '{}.SNP_Analysis.hgvs_fixed.txt').format(control_sample)
        cmd = ["munge", "hgvs_check", "--update", "--output-error-column",
               "-i", infile, "-o", outfile]
        subprocess.call(cmd)
        self.assertTrue(filecmp.cmp(outfile, expected))
