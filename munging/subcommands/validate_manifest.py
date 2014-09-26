"""
Usage:
validate_manifest
"""

import ConfigParser
import csv
import pandas
import re

def build_parser(parser):
        parser.add_argument('config',
        help='Config file')

def compare_manifest_to_pending(data):
    """Compare MRN, Accession, Name, Library, Well/Sample ID
    """
    pass

def parse_manifest(xls_file):
    """Gather the columns of interest:
    Run,well,barcode, capture, flowcell
    """
    xls = pandas.ExcelFile(xls_file)
    manifest = str([i for i in xls.sheet_names if re.search('manifest',i)]).strip('[]\'u')
    xls_data=xls.parse(manifest)#, index_col=0, na_values=['NA'])
    data = xls_data.to_dict()
    return data


#args = parser.parse_args()
def action(args):
    pending_list=pandas.read_csv(open(run_info['tracker'],'rU'))
    manifest_info=parse_manifest(run_info['manifest'])
    print "pending keys:", pending_list.keys()
    print "manifest keys:", manifest_info.keys()
    # import IPython
    # IPython.embed()
    config = ConfigParser.ConfigParser()
    config.read(args.config)
    run_info=config.defaults()



