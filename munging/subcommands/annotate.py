"""
Run annovar to generate standard set of annotations to bring into the DB using annotation_importer
"""
import sys
from collections import namedtuple
import os
#from os.path import join, isfile, dirname
import subprocess
import argparse
import csv
from munging.utils import munge_path, munge_pfx

def build_parser(parser):
    parser.add_argument('run_dir',
                        help='Directory where input file is located and where output files will be created')
    parser.add_argument('input_file', default=None,
                        help='Explicitly specify input file of variants in Annovar format')
    parser.add_argument('--library_dir', default='/mnt/disk2/com/Genomes/Annovar_files',
                        help='Directory holding Annovar library files')

    parser.add_argument('--annovar_bin', default='',
                        help='Location of the Annovar perl executables')

ANNOTATIONS = [('snp138',),  # dbsnp
               ('exac03',),  # ExAC 65000 exome allele frequency data 
               ('1000g2015aug_all',),  # 1000 genomes annotation:
               ('1000g2015aug_amr',),  # 1000 genomes (admmixed american) annotation:
               ('1000g2015aug_eur',),  # 1000 genomes (european) annotation:
               ('1000g2015aug_eas',),  # 1000 genomes (east asian) annotation:
               ('1000g2015aug_sas',),  # 1000 genomes (south asian) annotation:
               ('1000g2015aug_afr',),  # 1000 genomes (african) annotation:
               ('ljb26_all',),  # whole-exome SIFT, PolyPhen2 HDIV, PolyPhen2 HVAR, LRT, MutationTaster, MutationAssessor, FATHMM, MetaSVM, MetaLR, VEST, CADD, GERP++, PhyloP and SiPhy scores from dbNSFP version 2.6
               ('esp6500siv2_all',),  # alternative allele frequency in the NHLBI-ESP project with 6500 exomes, including the indel calls and the chrY calls. evs-all
               ('esp6500siv2_aa',),  # evs-african american
               ('esp6500siv2_ea',),  # evs-european
               ('cosmic70',),  # cosmic67
               ('clinvar_20150629',),  # CLINVAR database with Variant Clinical Significance (unknown, untested, non-pathogenic, probable-non-pathogenic, probable-pathogenic, pathogenic, drug-response, histocompatibility, other) and Variant disease name
               ('nci60',),  # NCI-60 human tumor cell line panel exome sequencing allele frequency data
               ('segdup', '--regionanno',),  # segdup annotation:               
               ('refGene', '--geneanno', ['--splicing_threshold','10', '--hgvs']),  # Gene level annotation:
               ]

# Named tuple for parsing annotation defs
AnnotInfo =  namedtuple('AnnotInfo', ['dbtype', 'anno_type', 'args'])
AnnotInfo.__new__.__defaults__ = ('-filter','')

def action(args):
    BUILDVER = 'hg19'
    ANNOVAR_VARIANTS = os.path.join(args.annovar_bin, 'annotate_variation.pl')
    GENERIC_DB = 'hg19_clinical_variants'

    pathinfo = munge_path(args.run_dir)
    internal_freq_file = '_'.join(['hg19',pathinfo['machine'],pathinfo['assay']])
    internal_cadd_file = '_'.join(['hg19','CADD',pathinfo['assay']])
    
    variants_file = args.input_file

    pfx_info = munge_pfx(os.path.basename(args.input_file))

    annots = [AnnotInfo(*a) for a in ANNOTATIONS]

    #Add the generics dbs to the annotation info
    if not pathinfo['assay'] == 'msi-plus':
        annots.append(AnnotInfo(dbtype='generic', anno_type='--genericdbfile', args=[internal_freq_file, '-filter']))
        GENERIC_DB='hg19_clinical_variants_msiplus'
    annots.append(AnnotInfo(dbtype='generic', anno_type='--genericdbfile', args=[internal_cadd_file, '-filter']))
    annots.append(AnnotInfo(dbtype='generic', anno_type='--genericdbfile', args=[GENERIC_DB, '-filter']))

    # run annotate_variants based on the spec in ANNOTATIONS
    for a in annots:
        annovar_cmd = [ANNOVAR_VARIANTS, 
                       '-dbtype', a.dbtype,
                       '--buildver', BUILDVER, 
                       a.anno_type ] \
                       + list(a.args) + \
                       ['--otherinfo', 
                        '--separate',
                        '-outfile', os.path.join(os.path.dirname(args.input_file), pfx_info['pfx']), 
                        variants_file, 
                        args.library_dir]
        cmd = filter(None, annovar_cmd)
        subprocess.check_call(cmd)
        if a.anno_type=='--genericdbfile':
            generic_file = os.path.join(os.path.dirname(args.input_file),pfx_info['pfx']+'.hg19_generic_dropped')
            gen_file_basename = a.args[0].replace(pathinfo['assay'],'').strip('_')
            if pathinfo['machine'] in a.args[0]:
                gen_file_basename = gen_file_basename.replace(pathinfo['machine']+'_'+pathinfo['assay'],'UW_freq')
            specific_file = os.path.join(os.path.dirname(args.input_file),pfx_info['pfx']+'.'+gen_file_basename+'_dropped')
            mvcmd=['mv' , generic_file, specific_file]
            subprocess.check_call(mvcmd)




