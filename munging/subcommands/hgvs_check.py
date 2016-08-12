"""
Use Mutaluzer.nl webservice to check HGVS c./p. variants for errors and optionally update them.

Usage:

munge hgvs_check --update input_file -o output_file

"""
from suds.client import Client

import argparse
import logging
import base64
import time
import csv
import sys
import re

log = logging.getLogger(__name__)

MUTALYZER_URL = 'https://mutalyzer.nl/services/?wsdl'


def build_parser(parser):
    parser.add_argument('-i', '--infile', type=argparse.FileType('r'), nargs='?',
                        default=sys.stdin, help='Input file')
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'), nargs='?',
                        default=sys.stdout, help='Output file')
    parser.add_argument('--update', action="store_true", default=False,
                        help='Use Mutalyzer to update HGVS',)
    parser.add_argument('--output-error-column', action="store_true", default=False,
                        help='Output additional column with any mutalyzer-reported'
                             'errors or warnings')
    parser.add_argument('--dna-column-name', type=str, default='c.',
                        help='Name of column with DNA/nucleotide HGVS annotation')
    parser.add_argument('--protein-column-name', type=str, default='p.',
                        help='Name of column with protein HGVS annotation')


def get_soap_connection(URL):
    c = Client(URL, cache=None)
    o = c.service
    return o


def query_mutalyzer(conn, batch_file_entries):
    encoded_data = base64.b64encode("\n".join(batch_file_entries) + "\n")
    log.info("Uploaded batch job data (%d rows)" % len(batch_file_entries))
    job_id = conn.submitBatchJob(data=encoded_data, process="NameChecker")
    time.sleep(5)
    while (conn.monitorBatchJob(job_id=job_id) > 0):
        log.info("Waiting for batch job with job_id %s to finish" % job_id)
        time.sleep(30)
    encoded_results = conn.getBatchJob(job_id=job_id)
    results = base64.b64decode(encoded_results).split("\n")
    return results


CODONS = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
          'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
          'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
          'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
          'TER': '*', 'XAA': 'X'}


def clean_protein_field(string):
    """
    Function to translate 3-letter IUPAC codes to 1-letter codes,
    as well as to remove parentheses from HGVS protein predictions.
    """
    if not string.startswith("p."):
        # not a protein annotation!
        return string
    string = re.sub(r'[()]', '', string)  # clean out parentheses
    if string == "p.=":
        # return blank if no protein change to be consistent with
        # ANNOVAR output
        return ""
    string = string[2:]
    outstring = "p."
    i = 0
    while i < len(string):
        if (i+3 < len(string)) and (string[i:i+3].upper() in CODONS):
            outstring += CODONS[string[i:i+3].upper()]
            i += 3
        else:
            outstring += str(string[i])
            i += 1
    return outstring


def action(args):
    delimiter = "\t"
    header = None
    batch_file_entries = []
    input_rows = []
    lookup_dict = {}  # dict input file mapping line number =>  mutalyzer input
    lookup_count = 0

    # read in data
    infile = csv.DictReader(args.infile, delimiter=delimiter)
    for line_num, line in enumerate(infile):
        # skip empty HGVS fields, but do some bookkeeping to associate
        # multiple HGVS per line. We therefore store a
        # dict mapping input line # to the position in the mutalyzer batch upload
        if line[args.dna_column_name] != '':
            hgvs_entries = line[args.dna_column_name].split(" ")
            lookup_dict[line_num] = range(lookup_count, lookup_count + len(hgvs_entries))
            lookup_count += len(hgvs_entries)
            batch_file_entries.extend(hgvs_entries)
        input_rows.append(line)

    # connect to mutalyzer and query batch file
    conn = get_soap_connection(MUTALYZER_URL)
    results = query_mutalyzer(conn, batch_file_entries)

    # split header from results; mutalyzer adds an extra blank row which we remove
    result_rows = results[1:-1]
    log.info("Received %d rows from mutalyzer web service" % (len(result_rows)))

    # ensure that Mutalyzer returns the same number of rows as we sent it
    if len(batch_file_entries) != len(result_rows):
        log.error("Wrong number of returned results! "
                  "Expecting %d, received %d" % (len(batch_file_entries), len(result_rows)))
        sys.exit(1)

    # Write output header
    if args.output_error_column:
        header = infile.fieldnames + ["mutalyzer_errors"]
    else:
        header = infile.fieldnames
    args.outfile.write(delimiter.join(header) + "\n")

    # Iterate through input rows and update/write lines to output
    for line_num, line in enumerate(input_rows):
        # check if there is a result
        ixs = lookup_dict.get(line_num, [])
        correct_hgvs_dna = []
        correct_hgvs_protein = []
        errors = []
        error_flag = False
        for ix in ixs:
            # get result line(s) and concat into error field
            result_row = result_rows[ix].rstrip().split("\t")
            if len(result_row) > 2:
                # construct new HGVS from mutalyzer output
                input_tx = result_row[0].split(":")[0]
                correct_hgvs_dna.append(input_tx + ":" + result_row[6])
                correct_hgvs_protein.append(clean_protein_field(result_row[7]))
            else:
                # mutalyzer has error, keep old hgvs
                error_flag = True
            errors.append(result_row[1])

        if args.update and not error_flag:
            old_hgvs_dna = line[args.dna_column_name]
            new_hgvs_dna = " ".join(correct_hgvs_dna)
            old_hgvs_protein = line[args.protein_column_name]
            new_hgvs_protein = " ".join(correct_hgvs_protein)
            if old_hgvs_dna != new_hgvs_dna:
                log.info("updating DNA HGVS: %s --> %s and protein HGVS: %s --> %s"
                         % (old_hgvs_dna, new_hgvs_dna, old_hgvs_protein, new_hgvs_protein))
                line[args.dna_column_name] = new_hgvs_dna
                line[args.protein_column_name] = new_hgvs_protein

        line["mutalyzer_errors"] = " ".join(errors)
        outrow = [line[field] for field in header]
        args.outfile.write(delimiter.join(outrow) + "\n")
