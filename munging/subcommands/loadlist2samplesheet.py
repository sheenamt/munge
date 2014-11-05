"""Converts Illumina SampleSheet CSV files to the run_info.yaml input file.
This allows running the analysis pipeline without Galaxy, using CSV input
files from Illumina SampleSheet or Genesifter.
"""
import os
import csv
import itertools
"""
Create SampleSheet for demux from LoadList provided by lab

Usage:

 munge samplesheet /path/to/load-list.csv

"""

# import argparse
import os
import sys
import subprocess
import IPython
import csv
from operator import itemgetter
from itertools import groupby

# parser = argparse.ArgumentParser()
def build_parser(parser):
    parser.add_argument('loadlist',
                    help="Absolute path to load list")

import difflib
import glob

# ## Create samplesheets

def _lane_detail_to_ss(fcid, ldetail, r):
    """Convert information about a lane into Illumina samplesheet output.
    """
    if len(ldetail['ControlWell'])>3:
        ldetail['SampleID']=('_').join([ldetail['PlateNumber'],ldetail['Well'],ldetail['Recipe'],ldetail['ControlWell']])
    else:
        ldetail['SampleID']=('_').join([ldetail['PlateNumber'],ldetail['Well'],ldetail['Recipe']])

    return [ldetail['FCID'], r, ldetail["SampleID"], 'hg19',
            ldetail["Index"], ldetail["Description"], "N", 
            ldetail["Recipe"], ldetail["Operator"],
            ldetail["SampleProject"]]
    
def write_sample_sheets(fcid, lane_details, out_dir=None):
    """Convert a flowcell into a samplesheet for demultiplexing.
    """
    fcid = fcid
    if out_dir is None:
        out_dir = run_folder
    out_file = os.path.join(out_dir, "%s.csv" % fcid)
    with open(out_file, "w") as out_handle:
        writer = csv.writer(out_handle)
        writer.writerow(["FCID", "Lane", "SampleID", "SampleRef", "Index",
                         "Description", "Control", "Recipe", "Operator", "SampleProject"])
        for ldetail in lane_details:
            for r in range(1,3):
                writer.writerow(_lane_detail_to_ss(fcid, ldetail, r))
    return out_file


def _get_flowcell_id(reader, require_single=True):
    """Retrieve the unique flowcell id represented in the SampleSheet.
    """
    fc_ids = set([x['FCID'] for x in reader])
    if require_single and len(fc_ids) > 1:
        raise ValueError("There are several FCIDs in the same samplesheet file: %s" % in_file)
    else:
        return fc_ids

def db_project_info(info):
    pass

def action(args):
    out_dir='test_output'
    reader=csv.DictReader(open(args.loadlist))
    lane_details = [row for row in reader]
    #SampleSheet.csv needs to be grouped by FCID
    write_sample_sheets(list(_get_flowcell_id(lane_details))[0], lane_details, out_dir)
    #Database needs info grouped by project
    db_project_info(lane_details)
