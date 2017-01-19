"""
Create SampleSheet for demux from LoadList provided by lab

Usage:

 munge loadlist2samplesheet /path/to/load-list.csv

"""

# import argparse
import os
import csv
import re

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

ASSAYS={'OPX':'OncoPlex',    
        'BRO':'ColoSeq',
        'EPI':'EpiPlex',
        'MRW':'MarrowSeq',
        'IMD':'ImmunoPlex',
        'CFF':'CellFreeFetal'}

def create_sample_project(ldetail):
    """Create sample project from Recipe and PlateNumber"""
    #Grab the assay based on the recipe, don't care about verion of assay for this part
    assay = [value for key,value in ASSAYS.items() if re.search(key, ldetail['Recipe'])][0]
    if ldetail['Description'].upper()=='KAPA':
        sample_project=assay+ldetail['Description'].upper()+ldetail['PlateNumber']
    elif ldetail['Description'].upper()=='STANDARD' or ldetail['Description'].upper()=='SURE SELECT':
        sample_project=assay+ldetail['PlateNumber']
    else:
        sample_project=assay+ldetail['Description']+ldetail['PlateNumber']
    # append assay version to the end - used in case of multiple assay versions
    sample_project +='-'+ldetail['Recipe']
    return sample_project

def _lane_detail_to_ss(fcid, ldetail, r):
    """Convert information about a lane into Illumina samplesheet output.
    FCID|PlateNumber|Well|Description|Control|Well|Recipe|Index
    """
    prefix=(ldetail['PlateNumber']+WELL_MAPPING[ldetail['Well']])

    if len(ldetail['ControlWell'])>4:
        ldetail['SampleID']=('-').join([prefix,ldetail['Well'],ldetail['Recipe'],ldetail['ControlWell']])
    else:
        ldetail['SampleID']=('-').join([prefix,ldetail['Well'],ldetail['Recipe']])

    ldetail["SampleProject"]=create_sample_project(ldetail)

    return [ldetail['FCID'], r, ldetail["SampleID"], 'hg19',
            ldetail["Index"], ldetail["Description"], "N", 
            ldetail["Recipe"],  ldetail["Operator"],
            ldetail["SampleProject"]]

def _lane_detail_to_signout(ldetail):
    """Convert information about a lane into Signout sheet
    SampleID |Accession|Patient Name| MRN
    """
    prefix=(ldetail['PlateNumber']+WELL_MAPPING[ldetail['Well']])

    if len(ldetail['ControlWell'])>5:
        ldetail['SampleID']=('_').join([prefix,ldetail['Well'],ldetail['Recipe'],ldetail['ControlWell']])
    else:
        ldetail['SampleID']=('_').join([prefix,ldetail['Well'],ldetail['Recipe']])

    ldetail["SampleProject"]=create_sample_project(ldetail)

    return ldetail["SampleID"], ldetail["Accession"],ldetail["Patient Name"],ldetail["MRN"]
    
def write_sample_sheet(fcid, lane_details, out_dir=None):
    """Convert a flowcell into a samplesheet for demultiplexing.
    """
    fcid = fcid
    out_file = open(os.path.join(out_dir, "%s.csv" % fcid), "w")
    signout=open(os.path.join(out_dir, "%s.signout.csv" % fcid),"w")
    writer = csv.writer(out_file)
    writer.writerow(["[Data]",])
    writer.writerow(["Sample_ID","Sample_Project","index"])
    so_writer = csv.writer(signout)
    so_writer.writerow(["SampleID", "Accession","Patient Name","MRN"])

    for ldetail in lane_details:
        so_writer.writerow(_lane_detail_to_signout(ldetail))
        info=_lane_detail_to_ss(fcid, ldetail, 1)
        writer.writerow([info[2],info[9],info[4]])
    return out_file

def check_control_vs_number_assays(lane_details):
    """Make sure there is NA12878 for each assay listed
    """
    controls = 0
    assays = set()
    for ldetail in lane_details:
        assays.add(ldetail['Recipe'])
        if ldetail['ControlWell']=='NA12878':
            controls = controls+1
    if len(assays) != controls:
        raise ValueError("There is not 1 NA12878 controls for each assay listed")

def _get_flowcell_id(reader, require_single=True):
    """Retrieve the unique flowcell id represented in the SampleSheet.
    """
    fc_ids = set([x['FCID'] for x in reader])
    if require_single and len(fc_ids) > 1:
        raise ValueError("There are several FCIDs in the same samplesheet file: %s" % fc_ids)
    else:
        return fc_ids

def action(args):
    out_dir='./'
    reader=csv.DictReader(open(args.loadlist,'rU'))
    #strip whitespace from header names in case tech used wrong template
    reader.fieldnames=[i.strip() for i in reader.fieldnames]
    #TGC wants to be sorted alphanumerically 
    lane_details = sorted(reader, key=lambda d: d['Well'])
    check_control_vs_number_assays(lane_details)
    #SampleSheet.csv needs to be grouped by FCID
    write_sample_sheet(list(_get_flowcell_id(lane_details))[0], lane_details, out_dir)
