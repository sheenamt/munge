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
import ssl


log = logging.getLogger(__name__)

EXTERNAL_MUTALYZER_URL = 'https://mutalyzer.nl/services/?wsdl'
INTERNAL_MUTALYZER_URL = 'https://mutalyzer.labmed.uw.edu/services/?wsdl'

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
    parser.add_argument('--server', type=str, default='redundant',
                        help='Specify internal, external, or redundant for server to use (redundant tries internal then external)')

def get_mutalyzer_connection(url):
    ''' Returns a service object if successful, otherwise returns a boolean value
        False.'''
    try:
        c = Client(url, cache=None)
        o = c.service
        return o
    except Exception as e:
        log.error(e)
        return False

def setup_mutalyzer_connection(server):
    ''' Contains the logic for which server to connect to.
        Returns a connection object if successful, a boolean
        False value otherwise. May need to test object type in calling 
        code. 
        Note: some clients do not recognize UW certs - thus
        ignoring SSL warnings when contacting internal UW server'''
    if server=='internal':
        ''' Establish internal server connection '''
        # Do not require SSL for internal server - UW certs sometimes not recognized
        ssl._create_default_https_context = ssl._create_unverified_context
        conn = get_mutalyzer_connection(INTERNAL_MUTALYZER_URL)
        if not isinstance(conn, bool):
            return conn
        else:
            log.error("Failed to connect to internal server, returning input file")
            return False

    elif server=='external':
        ''' Trying external server - only valid SSL connections are accepted.'''
        conn = get_mutalyzer_connection(EXTERNAL_MUTALYZER_URL)
        if not isinstance(conn, bool):
            log.info("Successful connection to external server established.")
            return conn
        else:
            log.error("Failed to connect to external server...returning input file.")
            return False

    else:
        ''' assume redundant - try internal, external, and internal ignore SSL if flagged.'''
        # Try to establish connections to both external and internal (ignoring internal SSL warning)
        external_conn = get_mutalyzer_connection(EXTERNAL_MUTALYZER_URL)

        # Create HTTPS context to ignore SSL warning when trying internal server
        ssl._create_default_https_context = ssl._create_unverified_context
        internal_conn = get_mutalyzer_connection(INTERNAL_MUTALYZER_URL)

        # If internal was successful use it, otherwise use external
        if not isinstance(internal_conn, bool):
            log.info("Successful connection to internal server established (Ignored SSL warning).")
            return internal_conn
        elif not isinstance(external_conn, bool):
            log.info("Successful connection to external server established.")
            return external_conn
        else:
            log.error("Internal and External connections failed...returning input file.")
            return False

def query_mutalyzer(conn, batch_file_entries):
    ''' Query mutalyzer'''
    # may need to evaluate efficacy of time calls here
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
        if (i+3 <= len(string)) and (string[i:i+3].upper() in CODONS):
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

    server = args.server
    conn = setup_mutalyzer_connection(server)

    # If we cannot establish a connection we simply return the original file
    if isinstance(conn, bool) and conn==False:
        log.info("Connection failed...returning input file with no amended mutalyzer values...")
        # Write input file back out
        for line in args.infile:
            args.outfile.write(line)
        sys.exit()

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
            if (old_hgvs_dna != new_hgvs_dna) or (old_hgvs_protein != new_hgvs_protein):
                log.info("updating DNA HGVS: %s --> %s and protein HGVS: %s --> %s"
                         % (old_hgvs_dna, new_hgvs_dna, old_hgvs_protein, new_hgvs_protein))
                line[args.dna_column_name] = new_hgvs_dna
                line[args.protein_column_name] = new_hgvs_protein

        line["mutalyzer_errors"] = " ".join(errors)
        outrow = [line[field] for field in header]
        args.outfile.write(delimiter.join(outrow) + "\n")

