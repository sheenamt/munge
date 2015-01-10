"""
Parsers for the top level summary files
pindel, cnv_exon, cnv_gene, snp, msi, quality, clin_flagged
"""
import os
import csv
import sys
import IPython

from itertools import count, groupby, chain, ifilter 
from operator import itemgetter

from munging import filters
from munging.utils import walker, munge_pfx
from _sqlite3 import Row

def parse_quality(files, specimens, annotation, prefixes, variant_keys):
    files = ifilter(filters.quality_analysis, files)
    files=sorted(files)    
    variant_keys = ['MEAN_TARGET_COVERAGE',]
 
    #sort the files so that the output in the workbook is sorted
    for pth in files:      
        pfx = munge_pfx(pth.fname)
        log_pfx=pfx['mini-pfx']
        prefixes.append(log_pfx)
        with open(os.path.join(pth.dir, pth.fname)) as fname:
            reader = csv.DictReader(fname, delimiter='\t')
            for row in reader:
                variant = tuple(k for k in variant_keys)
                specimens[variant][log_pfx] = row['MEAN_TARGET_COVERAGE']
                annotation[variant] = specimens[variant]

    fieldnames = variant_keys + prefixes
    return specimens, annotation, prefixes, fieldnames, variant_keys 

def parse_clin_flagged(files, specimens, annotation, prefixes, variant_keys):
    files = ifilter(filters.genotype_analysis, files)
    files=sorted(files)    
    variant_keys = ['Position','Ref_Base','Var_Base' ]
    
    #sort the files so that the output in the workbook is sorted
    for pth in files:
        pfx = munge_pfx(pth.fname)
        reads_pfx=pfx['mini-pfx']+'_Reads'
        prefixes.append(reads_pfx)
        with open(os.path.join(pth.dir, pth.fname)) as fname:
            reader = csv.DictReader(fname, delimiter='\t')
            for row in reader:
                variant = tuple(row[k] for k in variant_keys)
                specimens[variant][reads_pfx] = row['Reference_Reads']+'|'+row['Variant_Reads']
                annotation[variant] = row

    annotation_headers = [
        'Clinically_Flagged']
    fieldnames = variant_keys + annotation_headers + prefixes
    return specimens, annotation, prefixes, fieldnames, variant_keys            

def parse_msi(files, control_file, specimens, prefixes, variant_keys):
    files = ifilter(filters.msi_file_finder,files) 
    files=sorted(files)    
    control_info=csv.DictReader(control_file, delimiter='\t')
    control_info = sorted(control_info, key=itemgetter('Position'))
    control_info=[d for d in control_info]
    
    variant_keys = ['Position',]
    for pth in files:   
        pfx = munge_pfx(pth.fname)
        mini_pfx=pfx['mini-pfx']
        prefixes.append(mini_pfx)
        with open(os.path.join(pth.dir, pth.fname)) as fname:
            reader = csv.DictReader(fname, delimiter='\t')
            sample_msi = sorted(reader, key=itemgetter('Position'))
            for key, group in groupby(sample_msi, key=itemgetter('Position')):
                control_row=[d for d in control_info if d['Position']==key]
                variant = tuple(control_row[0][k] for k in variant_keys)    
                for sample_info in group:
                    if int(sample_info['Avg_read_depth']) >= 30:
                        value = float(control_row[0]['Ave']) + (2 * float(control_row[0]['Std']))
                        if int(sample_info['Number_Peaks']) >= value:
                            new_info = 1
                        else:
                            new_info = 0
                    else:           
                        new_info = None
                    specimens[variant][mini_pfx] = new_info

            for entry in specimens.items():
                for pfx in prefixes:
                    loci=0
                    if entry[1][pfx] is not None:
                       loci+=entry[1][pfx]
            specimens['Total_Loci'][pfx]=loci
    IPython.embed()
    fieldnames = variant_keys + list(prefixes)
    return specimens, prefixes, fieldnames, variant_keys            


def parse_pindel(files, specimens, annotation, prefixes, variant_keys):#SNP Specific    
    #Pindel Specific    
    files = ifilter(filters.pindel_analysis, files)
    files=sorted(files)    
    variant_keys = ['Position', 'Gene']
    annotation_headers = [
        'Gene_Region',
        'Event_Type',
        'Size',
        'Transcripts',
        ]
    #sort the files so that the output in the workbook is sorted
    files=sorted(files)
    for pth in files:
        pfx = munge_pfx(pth.fname)
        prefixes.append(pfx['mini-pfx'])
        with open(os.path.join(pth.dir, pth.fname)) as fname:
            reader = csv.DictReader(fname, delimiter='\t')
            for row in reader:
                variant = tuple(row[k] for k in variant_keys)
                specimens[variant][pfx['mini-pfx']] = row['Reads']
                annotation[variant] = row
                
    fieldnames = variant_keys + annotation_headers + prefixes
    return specimens, annotation, prefixes, fieldnames, variant_keys            

def parse_snp(files, specimens, annotation, prefixes, variant_keys):#SNP Specific   
    files = ifilter(filters.only_analysis, files)
    files = sorted(files)    
    variant_keys = ['Position', 'Ref_Base', 'Var_Base']
    annotation_headers = [
        'Gene',
        'Variant_Type',
        'Transcripts',
        'Clinically_Flagged',
        'Cosmic',
        'Segdup',
        'Polyphen',
        'Sift',
        'Mutation_Taster',
        'Gerp',
        'HiSeq_Freq',
        'HiSeq_Count',
        'MiSeq_Freq',
        'MiSeq_Count',
        '1000g_ALL',
        'EVS_esp6500_ALL',
        '1000g_AMR',
        'EVS_esp6500_AA',
        '1000g_EUR',
        'EVS_esp6500_EU',
        '1000g_ASN',
        '1000g_AFR']

    for pth in files:
        pfx = munge_pfx(pth.fname)
        print pth
        reads_pfx=pfx['mini-pfx']+'_Ref|Var'
        prefixes.append(reads_pfx)
        with open(os.path.join(pth.dir, pth.fname)) as fname:
            reader = csv.DictReader(fname, delimiter='\t')
            for row in reader:
                variant = tuple(row[k] for k in variant_keys)
                specimens[variant][reads_pfx] = row['Ref_Reads']+'|'+row['Var_Reads']
                annotation[variant] = row
    fieldnames = variant_keys + annotation_headers + prefixes
    return specimens, annotation, prefixes, fieldnames, variant_keys

def parse_cnv_exon(files, specimens, annotation, prefixes, variant_keys):
    
    files = ifilter(filters.cnv_exon_analysis, files)
    files = sorted(files)
    variant_keys = ['Position', 'Gene' ]
    #sort the files so that the output in the workbook is sorted
    for pth in files:
        print pth
        pfx = munge_pfx(pth.fname)
        log_pfx=pfx['mini-pfx']+'_Log'
        prefixes.append(log_pfx)
        with open(os.path.join(pth.dir, pth.fname)) as fname:
            reader = csv.DictReader(fname, delimiter='\t')
            for row in reader:
                variant = tuple(row[k] for k in variant_keys)
                specimens[variant][log_pfx] = row['Ave_Adjusted_Log_Ratio']
                annotation[variant] = row

    annotation_headers = [
        'Transcripts']
    
    fieldnames = variant_keys + annotation_headers + prefixes
    #print prefixes, fieldnames, variant_keys
    return specimens, annotation, prefixes, fieldnames, variant_keys

def parse_cnv_gene(files, specimens, annotation, prefixes, variant_keys):
    files = ifilter(filters.cnv_gene_analysis, files)
    files=sorted(files)
    variant_keys = ['Position', 'Gene' ]
    #sort the files so that the output in the workbook is sorted
    for pth in files:
        pfx = munge_pfx(pth.fname)
        log_pfx=pfx['mini-pfx']+'_Log'
        prefixes.append(log_pfx)
        with open(os.path.join(pth.dir, pth.fname)) as fname:
            reader = csv.DictReader(fname, delimiter='\t')
            for row in reader:
                variant = tuple(row[k] for k in variant_keys)
                specimens[variant][log_pfx] = row['Ave_Adjusted_Log_Ratio']
                annotation[variant] = row

    annotation_headers = [
        'Transcripts']
    fieldnames = variant_keys + annotation_headers + prefixes
    return specimens, annotation, prefixes, fieldnames, variant_keys
    
