"""
Create annotation of all variants in db (or only from GATK)

Usage:

munge db_annotation $SAVEPATH/variants.* -o $SAVEPATH/annotation_table.txt

"""
import re
import logging
import sys
import argparse
import csv
from collections import defaultdict
from os import path

from munging.annotation import get_location
from munging.subcommands.annovar_summary import map_headers

variant_headers = ['chr','start','stop','Ref_Base','Var_Base']

# (pattern, header_ids, var_key_ids)

file_types = {
#gatk files
    'variant_function': ({0: 'var_type_1',1: 'Gene',}, [2, 3, 4, 5, 6]),
    'exonic_variant_function':  ({1: 'var_type_2', 2:'Transcripts'}, [3, 4, 5, 6, 7]),
    'hg19_ALL.sites.2014_10_dropped': ({1: '1000g_ALL'}, [2, 3, 4, 5, 6]),
    'hg19_AMR.sites.2014_10_dropped': ({1: '1000g_AMR'}, [2, 3, 4, 5, 6]),
    'hg19_AFR.sites.2014_10_dropped': ({1: '1000g_AFR'}, [2, 3, 4, 5, 6]),
    'hg19_SAS.sites.2014_10_dropped': ({1: '1000g_SAS'}, [2, 3, 4, 5, 6]),
    'hg19_EAS.sites.2014_10_dropped': ({1: '1000g_EAS'}, [2, 3, 4, 5, 6]),
    'hg19_EUR.sites.2014_10_dropped': ({1: '1000g_EUR'}, [2, 3, 4, 5, 6]),
    'hg19_avsift_dropped': ({1: 'Sift'}, [2, 3, 4, 5, 6]),
    'hg19_cosmic70_dropped': ({1: 'Cosmic'}, [2, 3, 4, 5, 6]),
    'hg19_genomicSuperDups': ({0: 'Segdup'}, [2, 3, 4, 5, 6]),
    'hg19_ljb26_all_dropped': ({1: 'Polyphen'}, [2, 3, 4, 5, 6]),
    'hg19_ljb_gerp++_dropped': ({1: 'Gerp'}, [2, 3, 4, 5, 6]),
    'hg19_ljb_mt_dropped': ({1: 'Mutation_Taster'}, [2, 3, 4, 5, 6]),
    'hg19_esp6500si_all_dropped':({1: 'EVS_esp6500_ALL'}, [2, 3, 4, 5, 6]),
    'hg19_esp6500si_ea_dropped':({1: 'EVS_esp6500_EU'}, [2, 3, 4, 5, 6]),
    'hg19_esp6500si_aa_dropped':({1: 'EVS_esp6500_AA'}, [2, 3, 4, 5, 6]),
    'hg19_variants_dropped':({1:'Clinically_Flagged'}, [2, 3, 4, 5, 6]),
    'hg19_nci60_dropped':({1:'NCI60'},  [2, 3, 4, 5, 6]),
    'hg19_snp137_dropped':({1:'rsid_1'}, [2, 3, 4, 5, 6]),
    'hg19_clinvar_20140929_dropped': ({1: 'ClinVar'}, [2, 3, 4, 5, 6]),
}

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument(
        'infiles', action = 'append', nargs = '+',
        help='Input files')
    parser.add_argument(
        '-o','--outfile',
        help='Output file', default = sys.stdout,
        type = argparse.FileType('w'))
    parser.add_argument(
        '--strict', action = 'store_true', default = False,
        help = 'Exit with error if an input file has no match.')

    
def munge_gene_and_Transcripts(data):
    """
    Return modified values of (Gene, Transcripts). Note that
    this depends on 'Variant_Type' provided by munge_variant.
    """

    Transcripts = data.get('Transcripts')
    Gene = data.get('Gene', '')
    
    if not Gene or data['Variant_Type'] in ('upstream','downstream','intergenic','ncRNA_exonic'):
        Gene = ''
    elif '(' in Gene:        
        Gene, Gene_tail = data['Gene'].split('(', 1)
        if 'NM' in Gene_tail:
            # overrides value of Transcripts
            Transcripts = data['Gene']
                        
    return Gene, Transcripts

def action(args):

    (infiles, ) = args.infiles
    headers = ['Position'] + variant_headers[3:5] + [
        'Gene',
        'dbSNP_ID',
        'Variant_Type',
        'Transcripts',
        'Clinically_Flagged',
        'CADD',
        'ClinVar',
        'NCI60',
        'Cosmic',
        'Segdup',
        'Polyphen',
        'Sift',
        'Mutation_Taster',
        'Gerp',
        '1000g_ALL',
        'EVS_esp6500_ALL',
        '1000g_AMR',
        'EVS_esp6500_AA',
        '1000g_EUR',
        'EVS_esp6500_EU',
        '1000g_SAS',
        '1000g_EAS',
        '1000g_AFR']


    # accumulate data from all input files for each variant
    output = defaultdict(dict)
    for fname in infiles:
        _, file_type = path.basename(fname).split('.', 1)
        try:
            #if the file type matches one we want,
            #header ids are output columns
            #var_key_ids are chrm:str:stp:ref:var
            header_ids, var_key_ids = file_types[file_type]
        except KeyError:
            if re.search('dropped', fname):
                log.warning('no match: %s' % fname)
            if args.strict:
                sys.exit(1)
            continue

        for var_key, data in map_headers(fname, header_ids, var_key_ids):
            output[var_key].update(data)
                
    writer = csv.DictWriter(args.outfile,
                            fieldnames = headers,
                            quoting = csv.QUOTE_MINIMAL,
                            extrasaction = 'ignore',
                            delimiter='\t')

    writer.writeheader()
    sort_key = lambda row: [(row[k]) for k in ['chr','start','stop']]

    # write each row (with all data aggregated), modifying fields as necessary
    for data in sorted(output.values(), key = sort_key):
        # modify any specific fields here
        data['Variant_Type'] = data.get('var_type_2') if data.get('var_type_2','').strip() else data.get('var_type_1')
        data['Gene'], data['Transcripts'] = munge_gene_and_Transcripts(data)
        data['dbSNP_ID']=data.get('rsid_1') or data.get('rsid_2')
        data['1000g_ALL']=data.get('1000g_ALL') or -1
        data['1000g_AMR']=data.get('1000g_AMR') or -1
        data['1000g_SAS']=data.get('1000g_SAS') or -1
        data['1000g_EAS']=data.get('1000g_EAS') or -1
        data['1000g_AFR']=data.get('1000g_AFR') or -1
        data['1000g_EUR']=data.get('1000g_EUR') or -1
        data['EVS_esp6500_ALL']=data.get('EVS_esp6500_ALL') or -1
        data['EVS_esp6500_AA']=data.get('EVS_esp6500_AA') or -1
        data['EVS_esp6500_EU']=data.get('EVS_esp6500_EU') or -1
        writer.writerow(data)
