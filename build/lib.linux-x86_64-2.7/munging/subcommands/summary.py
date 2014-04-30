"""
Summarize output from Annovar and EVS

Usage:

   munge summary $SAVEPATH/$PFX.* -o $SAVEPATH/${PFX}_Analysis.txt;

"""

import itertools
import subprocess
import tempfile
import logging
import shutil
import os
import sys
import argparse
import re
import csv
from collections import defaultdict
from os import path

from munging.annotation import get_location, multi_split, split_string_in_two

variant_headers = ['chr','start','stop','Ref_Base','Var_Base']

# (pattern, header_ids, var_key_ids)

file_types = {
#gatk files

 'variant_function': ({0: 'var_type_1',
                            1: 'Gene',
                            7: 'Zygosity',
                            12: 'rsid_1',
                            8: 'GATK_Score'},
                           [2, 3, 4, 5, 6]),
 'exonic_variant_function':  ({1: 'var_type_2', 2:'Transcripts'}, [3, 4, 5, 6, 7]),
 'hg19_ALL.sites.2012_04_dropped': ({1: '1000g_ALL'}, [2, 3, 4, 5, 6]),
 'hg19_AMR.sites.2012_04_dropped': ({1: '1000g_AMR'}, [2, 3, 4, 5, 6]),
 'hg19_AFR.sites.2012_04_dropped': ({1: '1000g_AFR'}, [2, 3, 4, 5, 6]),
 'hg19_ASN.sites.2012_04_dropped': ({1: '1000g_ASN'}, [2, 3, 4, 5, 6]),
 'hg19_EUR.sites.2012_04_dropped': ({1: '1000g_EUR'}, [2, 3, 4, 5, 6]),
 'hg19_avsift_dropped': ({1: 'Sift'}, [2, 3, 4, 5, 6]),
 'hg19_cosmic65': ({1: 'Cosmic'}, [2, 3, 4, 5, 6]),
 'hg19_genomicSuperDups': ({0: 'Segdup'}, [2, 3, 4, 5, 6]),
 'hg19_ljb_all_dropped': ({1: 'Polyphen'}, [2, 3, 4, 5, 6]),
 'hg19_ljb_gerp++_dropped': ({1: 'Gerp'}, [2, 3, 4, 5, 6]),
 'hg19_ljb_mt_dropped': ({1: 'Mutation_Taster'}, [2, 3, 4, 5, 6]),
 'hg19_esp6500si_all_dropped':({1: 'EVS_esp6500_ALL'}, [2, 3, 4, 5, 6]),
 'hg19_esp6500si_ea_dropped':({1: 'EVS_esp6500_EU'}, [2, 3, 4, 5, 6]),
 'hg19_esp6500si_aa_dropped':({1: 'EVS_esp6500_AA'}, [2, 3, 4, 5, 6]),
 'hg19_miseq_dropped':({1: 'Mi_Freq_list'}, [2, 3, 4, 5, 6]),
 'hg19_hiseq_dropped':({1: 'Hi_Freq_list'}, [2, 3, 4, 5, 6]),
 'hg19_variants_dropped':({1:'Clinically_Flagged'}, [2, 3, 4, 5, 6]),
 'hg19_nci60_dropped':({1:'NCI60'},  [2, 3, 4, 5, 6]),

# varscanINDELS files
  'varscan.variant_function': ({0: 'var_type_1',
                                1: 'Gene',
                                7: 'Reads',
                                8:'Variant_Phred'},
                                [2, 3, 4, 5, 6]),
 'varscan.exonic_variant_function': ({1: 'var_type_2', 2:'Transcripts'}, [3, 4, 5, 6, 7]),
 'varscan.hg19_AFR.sites.2012_04_dropped': ({1: '1000g_AFR'}, [2, 3, 4, 5, 6]),
 'varscan.hg19_ALL.sites.2012_04_dropped': ({1: '1000g_ALL'}, [2, 3, 4, 5, 6]),
 'varscan.hg19_AMR.sites.2012_04_dropped': ({1: '1000g_AMR'}, [2, 3, 4, 5, 6]),
 'varscan.hg19_ASN.sites.2012_04_dropped': ({1: '1000g_ASN'}, [2, 3, 4, 5, 6]),
 'varscan.hg19_EUR.sites.2012_04_dropped': ({1: '1000g_EUR'}, [2, 3, 4, 5, 6]),
 'varscan.hg19_avsift_dropped':({1: 'Sift'}, [2, 3, 4, 5, 6]),
 'varscan.hg19_cosmic65_dropped':({1: 'Cosmic'}, [2, 3, 4, 5, 6]),
 'varscan.hg19_ljb_all_dropped': ({1: 'Poylphen'}, [2, 3, 4, 5, 6]),
 'varscan.hg19_ljb_gerp++_dropped': ({1: 'Gerp'}, [2, 3, 4, 5, 6]),
 'varscan.hg19_ljb_mt_dropped': ({1: 'Mutation_Taster'}, [2, 3, 4, 5, 6]),
 'varscan.hg19_snp137_dropped':({1: 'rsid_2'}, [2, 3, 4, 5, 6]),
 'varscan.hg19_esp6500si_all_dropped':({1: 'EVS_esp6500_ALL'}, [2, 3, 4, 5, 6]),
 'varscan.hg19_esp6500si_ea_dropped':({1: 'EVS_esp6500_EU'}, [2, 3, 4, 5, 6]),
 'varscan.hg19_esp6500si_aa_dropped':({1: 'EVS_esp6500_AA'}, [2, 3, 4, 5, 6]),
 'varscan.hg19_genomicSuperDups':({0: 'Segdup'}, [2, 3, 4, 5, 6]),
 'varscan.hg19_miseq_dropped':({1: 'Mi_Freq_list'}, [2, 3, 4, 5, 6]),
 'varscan.hg19_hiseq_dropped':({1: 'Hi_Freq_list'}, [2, 3, 4, 5, 6]),
 'varscan.hg19_variants_dropped':({1: 'Clinically_Flagged'}, [2, 3, 4, 5, 6]),
 'varscan.hg19_nci60_dropped':({1:'NCI60'},  [2, 3, 4, 5, 6]),

 #varscanSNP files
 'varscanSNP.variant_function': ({0: 'var_type_1',
                                 1: 'Gene',
                                  7: 'Reads',
                                  8:'Variant_Phred'},
                                [2, 3, 4, 5, 6]),
 'varscanSNP.exonic_variant_function': ({1: 'var_type_2', 2:'Transcripts'}, [3, 4, 5, 6, 7]),
 'varscanSNP.hg19_AFR.sites.2012_04_dropped': ({1: '1000g_AFR'}, [2, 3, 4, 5, 6]),
 'varscanSNP.hg19_ALL.sites.2012_04_dropped': ({1: '1000g_ALL'}, [2, 3, 4, 5, 6]),
 'varscanSNP.hg19_AMR.sites.2012_04_dropped': ({1: '1000g_AMR'}, [2, 3, 4, 5, 6]),
 'varscanSNP.hg19_ASN.sites.2012_04_dropped': ({1: '1000g_ASN'}, [2, 3, 4, 5, 6]),
 'varscanSNP.hg19_EUR.sites.2012_04_dropped': ({1: '1000g_EUR'}, [2, 3, 4, 5, 6]),
 'varscanSNP.hg19_avsift_dropped':({1: 'Sift'}, [2, 3, 4, 5, 6]),
 'varscanSNP.hg19_cosmic65_dropped':({1: 'Cosmic'}, [2, 3, 4, 5, 6]),
 'varscanSNP.hg19_ljb_all_dropped':({1: 'Polyphen'}, [2, 3, 4, 5, 6]),
 'varscanSNP.hg19_ljb_gerp++_dropped': ({1: 'Gerp'}, [2, 3, 4, 5, 6]),
 'varscanSNP.hg19_ljb_mt_dropped':({1:'Mutation_Taster'}, [2, 3, 4, 5, 6]),
 'varscanSNP.hg19_snp137_dropped':({1:'rsid_2'}, [2, 3, 4, 5, 6]),
 'varscanSNP.hg19_esp6500si_all_dropped':({1:'EVS_esp6500_ALL'}, [2, 3, 4, 5, 6]),
 'varscanSNP.hg19_esp6500si_ea_dropped':({1:'EVS_esp6500_EU'}, [2, 3, 4, 5, 6]),
 'varscanSNP.hg19_esp6500si_aa_dropped':({1:'EVS_esp6500_AA'}, [2, 3, 4, 5, 6]),
 'varscanSNP.hg19_genomicSuperDups':({0:'Segdup'}, [2, 3, 4, 5, 6]),
 'varscanSNP.hg19_miseq_dropped':({1:'Mi_Freq_list'}, [2, 3, 4, 5, 6]),
 'varscanSNP.hg19_hiseq_dropped':({1:'Hi_Freq_list'}, [2, 3, 4, 5, 6]),
 'varscanSNP.hg19_variants_dropped':({1:'Clinically_Flagged'}, [2, 3, 4, 5, 6]),
 'varscanSNP.hg19_nci60_dropped':({1:'NCI60'},  [2, 3, 4, 5, 6]),
}



log = logging.getLogger(__name__)

def munge_gene_and_Transcripts(data, preferred_trans):
    """
    Return modified values of (Gene, Transcripts). Note that
    this depends on 'Variant_Type' provided by munge_variant.
    """
    c,p="",""
    Transcripts = data.get('Transcripts')
    c,p=munge_transcript(Transcripts, preferred_trans)
    Gene = data.get('Gene', '')
    if not Gene or data['Variant_Type'] in ('upstream','downstream','intergenic','ncRNA_exonic'):
        Gene = ''
    elif '(' in Gene:        
        Gene, Gene_tail = data['Gene'].split('(', 1)
        if 'NM' in Gene_tail:
            # overrides value of Transcripts
            Transcripts = data['Gene']
            c,p=munge_transcript(Transcripts, preferred_trans)
    return Gene, Transcripts, c, p
        
def munge_transcript(data, preferred_trans):
    """
    Return HGVS correct transcript annotations
    Filtered with a preferred transcript list
    NM_006772.1:c.1713G>A
    """
    c,p="",""
    if not data==None:
        #Split incoming trans, strip the trailing )
        data=multi_split(data.strip(')'), ',(')
        for d in data:
            #Split the actual transcript info which is colon separated
            x=d.split(':')
            if len(x)==5:
                #x=['MTHFR', 'NM_005957', 'exon8', 'c.1305C>T', 'p.F435F']
                for pr in preferred_trans:
                    if re.search(x[1]+".", pr):
                        c=c+pr+':'+x[3]+' '
                        p=p+pr+':'+x[4]+' '
            elif len(x)==3:
                #x=['NM_001204171', 'exon5', 'c.343+9C>T']
                for pr in preferred_trans:
                    if re.search(x[0]+".", pr):
                        c=c+pr+':'+x[2]+' '
    return c.strip(), p.strip()

def map_headers(fname, header_ids, variant_idx):
    """
    Gets header(s) and info from each file
    """
    # dict mapping variant keys to columns
    with open(fname, 'rU') as infile:
        reader = csv.reader(infile, delimiter='\t')
        for row in reader:

            #adds the value from each position (2,3,4,5,6) to the variant_id
            variant_id = tuple(row[i] for i in variant_idx)

            #create dict of variants {2:'chr', 3:'start', 4:'stop, 5:'ref', 6:'var'}
            variant_dict = dict(zip(variant_headers, variant_id))
            
            #create dict of data from this file; key is included only
            #if value is not the empty string
            data = dict((key, row[i]) for i, key in header_ids.items() if row[i].strip())
            data['Position'] = get_location(**variant_dict)
            
            # if 'var_type_2' in header_ids.values():
            #     print variant_id, fname, data['var_type_2']
            
            yield (variant_id, dict(data, **variant_dict))

def get_allele_freq(data):
    """
    Return allele frequency of var_reads/ref_reads
    """
    if int(data['Var_Reads'])==-1 and int(data['Ref_Reads'])==-1:
        freq='NA'
    else:    
        total=int(data['Var_Reads'])+int(data['Ref_Reads'])
        freq=int(data['Var_Reads'])/float(total)
        freq="{0:.2f}".format(freq)
    return freq

def build_parser(parser):
    parser.add_argument('captured_genes', type=argparse.FileType('rU'), 
        help= 'Capture genes file with RefSeq in second column')
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

    
def action(args):

    (infiles, ) = args.infiles
    cap_file=args.captured_genes
    cap_data=csv.DictReader(cap_file, delimiter='\t')
    preferred_trans=[]
    annovar_trans=""

    for row in cap_data:
        for i in row['RefSeq'].split('/'):
            preferred_trans.append(i)

    headers = ['Position'] + variant_headers[3:5] + [
        'Clinically_Flagged',
        'Variant_Type',
        'HiSeq_Freq',
        '1000g_ALL',
        'EVS_esp6500_ALL',
        'Gene',
        'p.',
        'c.',
        'Faves_Y/N',
        'Ref_Reads',
        'Var_Reads',
        'Allele_Freq',
        'Variant_Phred',
        'Cosmic',
        'Polyphen',
        'Sift',
        'Mutation_Taster',
        'Gerp',
        '1000g_AMR',
        '1000g_EUR',
        '1000g_ASN',
        '1000g_AFR',
        'EVS_esp6500_AA',
        'EVS_esp6500_EU',
        'Transcripts',
        'Zygosity',
        'Segdup',
        'NCI60',
        'dbSNP_ID',
        'HiSeq_Count',
        'MiSeq_Freq',
        'MiSeq_Count',
        'GATK_Score'
        ]


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
        # # modify any specific fields here
        data['Variant_Type'] = data.get('var_type_2') if data.get('var_type_2','').strip() else data.get('var_type_1')
        data['Gene'], data['Transcripts'],data['c.'], data['p.'] = munge_gene_and_Transcripts(data, preferred_trans)
        data['dbSNP_ID']=data.get('rsid_1') or data.get('rsid_2')
        data['Variant_Phred']=data.get('Variant_Phred') 
        data['1000g_ALL']=data.get('1000g_ALL') or -1
        data['1000g_AMR']=data.get('1000g_AMR') or -1
        data['1000g_ASN']=data.get('1000g_ASN') or -1
        data['1000g_AFR']=data.get('1000g_AFR') or -1
        data['1000g_EUR']=data.get('1000g_EUR') or -1
        data['EVS_esp6500_ALL']=data.get('EVS_esp6500_ALL') or -1
        data['EVS_esp6500_AA']=data.get('EVS_esp6500_AA') or -1 
        data['EVS_esp6500_EU']=data.get('EVS_esp6500_EU') or -1 
        data['Ref_Reads'], data['Var_Reads']=split_string_in_two(data.get('Reads'))
        data['MiSeq_Freq'], data['MiSeq_Count']=split_string_in_two(data.get('Mi_Freq_list'))
        data['HiSeq_Freq'], data['HiSeq_Count']=split_string_in_two(data.get('Hi_Freq_list'))
        data['Allele_Freq']=get_allele_freq(data)
        writer.writerow(data)
