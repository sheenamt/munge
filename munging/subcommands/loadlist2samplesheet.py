"""
Create SampleSheet for demux from LoadList provided by lab

Usage:

 munge loadlist2samplesheet /path/to/load-list.csv

"""

# import argparse
import os
import sys
import subprocess
import csv
from operator import itemgetter
from itertools import groupby

# parser = argparse.ArgumentParser()
def build_parser(parser):
    parser.add_argument('loadlist',
                    help="Absolute path to load list")

WELL_MAPPING={'A01':'01',
              'B01':'02',
              'C01':'03',
              'D01':'04',
              'E01':'05',
              'F01':'06',
              'G01':'07',
              'H01':'08',
              'A02':'09',
              'B02':'10',
              'C02':'11',
              'D02':'12',
              'E02':'13',
              'F02':'14',
              'G02':'15',
              'H02':'16',
              'A03':'17',
              'B03':'18',
              'C03':'19',
              'D03':'20',
              'E03':'21',
              'F03':'22',
              'G03':'23',
              'H03':'24',
              'A04':'25',
              'B04':'26',
              'C04':'27',
              'D04':'28',
              'E04':'29',
              'F04':'30',
              'G04':'31',
              'H04':'32',
              'A05':'33',
              'B05':'34',
              'C05':'35',
              'D05':'36',
              'E05':'37',
              'F05':'38',
              'G05':'39',
              'H05':'40',
              'A06':'41',
              'B06':'42',
              'C06':'43',
              'D06':'44',
              'E06':'45',
              'F06':'46',
              'G06':'47',
              'H06':'48',
              'A07':'49',
              'B07':'50',
              'C07':'51',
              'D07':'52',
              'E07':'53',
              'F07':'54',
              'G07':'55',
              'H07':'56',
              'A08':'57',
              'B08':'58',
              'C08':'59',
              'D08':'60',
              'E08':'61',
              'F08':'62',
              'G08':'63',
              'H08':'64',
              'A09':'65',
              'B09':'66',
              'C09':'67',
              'D09':'68',
              'E09':'69',
              'F09':'70',
              'G09':'71',
              'H09':'72',
              'A10':'73',
              'B10':'74',
              'C10':'75',
              'D10':'76',
              'E10':'77',
              'F10':'78',
              'G10':'79',
              'H10':'80',
              'A11':'81',
              'B11':'82',
              'C11':'83',
              'D11':'84',
              'E11':'85',
              'F11':'86',
              'G11':'87',
              'H11':'88',
              'A12':'89',
              'B12':'90',
              'C12':'91',
              'D12':'92',
              'E12':'93',
              'F12':'94',
              'G12':'95',
              'H12':'96'}

ASSAYS={'OPXv4':'OncoPlex',    
        'OPXv3':'OncoPlex',    
        'BROv7':'ColoSeq',
        'BROv8':'ColoSeq',
        'BROv6':'ColoSeq',
        'EPIv1':'EpiPlex',
        'MRWv3':'MarrowSeq',
        'IMMv1':'ImmunoPlex'}

def create_sample_project(ldetail):
    """Create sample project from Recipe and PlateNumber"""
    if ldetail['Description'].upper() == 'KAPA':
        sample_project=ASSAYS[ldetail['Recipe']]+ldetail['Description'].upper()+ldetail['PlateNumber']
    else:
        sample_project=ASSAYS[ldetail['Recipe']]+ldetail['PlateNumber']
    return sample_project

def _lane_detail_to_ss(fcid, ldetail, r):
    """Convert information about a lane into Illumina samplesheet output.
    FCID|PlateNumber|Well|Description|Control|Well|Recipe|Index
    """
    prefix=(ldetail['PlateNumber']+WELL_MAPPING[ldetail['Well']])

    if len(ldetail['ControlWell'])>3:
        ldetail['SampleID']=('_').join([prefix,ldetail['Well'],ldetail['Recipe'],ldetail['ControlWell']])
    else:
        ldetail['SampleID']=('_').join([prefix,ldetail['Well'],ldetail['Recipe']])

    ldetail["SampleProject"]=create_sample_project(ldetail)

    return [ldetail['FCID'], r, ldetail["SampleID"], 'hg19',
            ldetail["Index"], ldetail["Description"], "N", 
            ldetail["Recipe"],  ldetail["Operator"],
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
    out_dir='./'
    reader=csv.DictReader(open(args.loadlist))
    lane_details = [row for row in reader]
    
    #SampleSheet.csv needs to be grouped by FCID
    write_sample_sheets(list(_get_flowcell_id(lane_details))[0], lane_details, out_dir)
    #Database needs info grouped by project
    db_project_info(lane_details)
