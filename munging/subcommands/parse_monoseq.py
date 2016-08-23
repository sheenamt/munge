"""
Run MonoSeq and parse output into readable format

"""

import subprocess
import tempfile
import logging
import shutil
import os
import xml.etree.ElementTree as ET
import csv

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('infile', 
        help='A required input file')
    parser.add_argument('xmlfile', 
        help='A required xml file that describes homopolyer site')
    parser.add_argument('outfile', 
        help='A required outputfile')
def action(args):

    tree = ET.parse(args.xmlfile)
    root = tree.getroot()
    homopoly = root[0][1].text[0]
    
    with open(args.infile) as f:
        lines=f.read().splitlines()
        polyt_info = lines[1].split(' ')
        #Monoseq also prints info for 0 length,so add one to this value
        num_polys = int(polyt_info[1])+1
    
    #header
    call_freqs = []
    output = {}
    for i in range(0,num_polys):
        call_freqs.append(str(i)+homopoly)

    call_info = dict(zip(call_freqs, polyt_info[-num_polys:]))
    output['COVERAGE']=polyt_info[0]
    output['SITE']=args.xmlfile.strip('data/ .xml')

    for key in sorted(call_info.iterkeys(), reverse=True):
        if float(call_info[key])>0.0:
            output[key]=call_info[key]

    print output
    header = sorted(output.keys(), reverse=True)
    print 'header:', header
#    writer = csv.DictWriter(open(args.outfile, 'w'),
    writer = csv.DictWriter(open(args.outfile, 'w'),
                        fieldnames = header, 
                        quoting=csv.QUOTE_NONE,
                        delimiter='\t')

    writer.writeheader()
    writer.writerow(output)




