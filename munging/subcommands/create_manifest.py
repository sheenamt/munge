"""
Usage:
create_manifest
"""

import ConfigParser
import csv

config = ConfigParser.ConfigParser()
config.read('config')
run_info=config.defaults()

def build_parser(parser):
    parser.add_argument('path',
                        help='Path to analysis files')

def parse_tracker(input_file):
    """Given tracker CSV, return relevant dict"""
    pending = csv.DictReader(open(input_file, 'rU'))

def create_manifest_dict(run_info):
    """Its so specific, it needs its own function"""
    pass

def write_manifest():
    """Create the excel manifest"""
    pass

def write_manifest():
    """Create the excel manifest"""
    pass

def action(args):
    pass
#args = parser.parse_args()
# pending_list=parse_tracker(run_info['input'])
# manifest_info=create_manifest_dict(run_info)



