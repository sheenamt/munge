"""
Parsers for the top level summary files
pindel, cnv_exon, cnv_gene, snp, msi, quality, clin_flagged
"""
import os
import csv
import sys
import copy
import re

from itertools import count, groupby, chain, ifilter , izip_longest
from operator import itemgetter

from munging import filters

"""Each function parses a group of sample files for desired information,
grouping based on the variant_keys list,
some include additional annotation headers,
sample counts, and scores calculated based on counts
"""

def shorten_name(fname):
    """ Shorten file name if needed """
    pfx=re.split('_|-', fname.split('.')[0])
    if 'NA12878' in pfx:
        pfx = '_'.join([pfx[0],pfx[3]])
        return pfx
    elif len(pfx) >2:
        pfx = '_'.join(pfx[:2])
        return pfx
    else:
        return fname.split('.')[0]

def parse_quality(files, specimens, annotation, prefixes, variant_keys):
    """ Parse the sample quality analysis file, from hs_metrics"""
    files = ifilter(filters.quality_analysis, files)
    files=sorted(files)    
    variant_keys = ['MEAN_TARGET_COVERAGE']
 
    #sort the files so that the output in the workbook is sorted
    for pth in files:      
        pfx = shorten_name(pth.fname)
        prefixes.append(pfx)
        with open(os.path.join(pth.dir, pth.fname)) as fname:
            reader = csv.DictReader(fname, delimiter='\t')
            for row in reader:
                variant = tuple(k for k in variant_keys)
                specimens[variant][pfx] = row['MEAN_TARGET_COVERAGE']
                annotation[variant] = specimens[variant]

    fieldnames = variant_keys + prefixes
    return specimens, annotation, prefixes, fieldnames, variant_keys 

def parse_clin_flagged(files, specimens, annotation, prefixes, variant_keys):
    """Parse the Genotype output, which is the reads of clin_flagged found"""
    files = ifilter(filters.genotype_analysis, files)
    files=sorted(files)    
    variant_keys = ['Position','Ref_Base','Var_Base' ]
    #sort the files so that the output in the workbook is sorted
    for pth in files:
        pfx = shorten_name(pth.fname)
        reads_pfx=pfx+'_Variants'
        prefixes.append(reads_pfx)
        with open(os.path.join(pth.dir, pth.fname)) as fname:
            reader = csv.DictReader(fname, delimiter='\t')
            for row in reader:
                variant = tuple(row[k] for k in variant_keys)
                specimens[variant][reads_pfx]=row['Variant_Reads']
                annotation[variant] = row

    annotation_headers = ['Clinically_Flagged']
    fieldnames = variant_keys + annotation_headers + prefixes
    return specimens, annotation, prefixes, fieldnames, variant_keys            

def parse_msi_flagged(files, specimens, annotation, prefixes, variant_keys):
    """Parse the Genotype output, which is the reads of clin_flagged found"""
    files = ifilter(filters.genotype_analysis, files)
    files=sorted(files)    
    variant_keys = ['Position','Ref_Base','Var_Base' ]
    #sort the files so that the output in the workbook is sorted
    for pth in files:
        pfx = shorten_name(pth.fname)
        reads_pfx=pfx+'_Variants|Total'
        status_pfx=pfx+'_Status'
        prefixes.append(reads_pfx)
        prefixes.append(status_pfx)

        with open(os.path.join(pth.dir, pth.fname)) as fname:
            reader = csv.DictReader(fname, delimiter='\t')
            for row in reader:
                variant = tuple(row[k] for k in variant_keys)
                try:
                    frac = "{0:.4f}".format(float(row['Variant_Reads'])/float(row['Valid_Reads']))
                except ZeroDivisionError:
                    frac = '0'
                specimens[variant][reads_pfx]=row['Variant_Reads']+'|'+row['Valid_Reads']
                if int(row['Valid_Reads']) >= 100:
                    if float(frac) >= 0.02:
                        specimens[variant][status_pfx]='POS'
                    elif float(frac) <= 0.01:
                        specimens[variant][status_pfx]='NEG'
                    elif 0.01 < float(frac) < 0.02 :
                        specimens[variant][status_pfx]='IND'
                else:
                    specimens[variant][status_pfx]='REVIEW'
                annotation[variant] = row
    annotation_headers = ['Clinically_Flagged']
    fieldnames = variant_keys + annotation_headers + prefixes
    return specimens, annotation, prefixes, fieldnames, variant_keys            



def parse_pindel(files, specimens, annotation, prefixes, variant_keys):
    """Parse the pindel analysis file, give total counts of samples with site"""
    #Grab just the pindel files
    files = ifilter(filters.pindel_analysis, files)
    #sort the files so that the output in the workbook is sorted
    files=sorted(files)    
    #List of keys to group samples by
    variant_keys = ['Position', 'Gene']
    #Other annotation to keep 
    annotation_headers = [
        'Gene_Region',
        'Event_Type',
        'Size',
        'Transcripts'
        ]

    #Go through all the files
    for pth in files:
        pfx = shorten_name(pth.fname)
        #Concatenate the pfx to human readable
        prefixes.append(pfx)
        with open(os.path.join(pth.dir, pth.fname)) as fname:
            reader = csv.DictReader(fname, delimiter='\t')
            for row in reader:
                variant = tuple(row[k] for k in variant_keys)
                #Update the specimen dict for this variant, for this pfx, report the Reads found
                specimens[variant][pfx] = row['Reads']
                annotation[variant] = row

    #Update the specimen dict for this variant, count samples present
    for key, value in specimens.iteritems():
        specimens[key]['Count']=len(value)

    #Add 'Count' to prefixes for correct dict zipping/printing    
    prefixes.append('Count')
    fieldnames = variant_keys + annotation_headers + prefixes
    return specimens, annotation, prefixes, fieldnames, variant_keys            

def parse_snp(files, specimens, annotation, prefixes, variant_keys):#SNP Specific   
    """Parse the snp output file, give ref|var read counts per sample"""
    files = ifilter(filters.snp_analysis, files)
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
        'NextSeq_Freq',
        'NextSeq_Count',
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
        pfx = shorten_name(pth.fname)
        reads_pfx=pfx+'_Ref|Var'
        prefixes.append(reads_pfx)
        with open(os.path.join(pth.dir, pth.fname)) as fname:
            reader = csv.DictReader(fname, delimiter='\t')
            for row in reader:
                variant = tuple(row[k] for k in variant_keys)
                specimens[variant][reads_pfx] = row['Ref_Reads']+'|'+row['Var_Reads']
                annotation[variant] = row

    #Update the specimen dict for this variant, count samples present
    for key, value in specimens.iteritems():
        specimens[key]['Count']=len(value)

    #Add 'Count' to prefixes for correct dict zipping/printing    
    prefixes.append('Count')
    fieldnames = variant_keys + annotation_headers + prefixes
    return specimens, annotation, prefixes, fieldnames, variant_keys            


def parse_cnv_exon(files, specimens, annotation, prefixes, variant_keys):
    """Parse the cnv_exon output, give ave_log_ratio"""
    files = ifilter(filters.cnv_exon_analysis, files)
    files = sorted(files)
    variant_keys = ['Position', 'Gene' ]
    #sort the files so that the output in the workbook is sorted
    for pth in files:
        pfx = shorten_name(pth.fname)
        log_pfx=pfx+'_Log'
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
    """Parse the cnv_genes output, give ave_log_ratio"""
    files = ifilter(filters.cnv_gene_analysis, files)
    files=sorted(files)
    variant_keys = ['Position', 'Gene' ]
    #sort the files so that the output in the workbook is sorted
    for pth in files:
        pfx = shorten_name(pth.fname)
        log_pfx=pfx+'_Log'
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
    

def parse_hsmetrics(lines, variant_keys):
    """
    Create human readable output from picard hsmetrics file
    """

    variant_keys=['MEAN_TARGET_COVERAGE',
                  'PF_READS',
                  'PF_UNIQUE_READS',
                  'PCT_PF_UQ_READS',
                  'PF_UQ_READS_ALIGNED',
                  'PCT_SELECTED_BASES',
                  'PCT_OFF_BAIT',
                  'PCT_USABLE_BASES_ON_TARGET',
                  'ZERO_CVG_TARGETS_PCT',
                  'AT_DROPOUT',
                  'GC_DROPOUT']

    #filter out lines that start with # or are just new lines    
    datalines=filter(lambda x: '#' not in x, lines)
    datalines=filter(lambda x: x!='\n' , datalines)

    #then strip newlines from the lines that contain data
    datalines=map(lambda x: x.strip() , datalines)

    keys=datalines[0].split('\t')
    values=datalines[1].split('\t')

    #metrics_dict=dict(zip(keys,values))
    metrics_dict = dict(izip_longest(keys,values, fillvalue='NA'))

    #Only print the keys we care about:
    output_dict = dict(zip(variant_keys,itemgetter(*variant_keys)(metrics_dict)))
    return output_dict, variant_keys

def parse_qualitymetrics(lines, variant_keys):
    """
    Create human readable output from picard hsmetrics file
    """

    variant_keys=['PERCENT_DUPLICATION',     
                  'READ_PAIRS_EXAMINED',
                  'READ_PAIR_DUPLICATES',
                  'UNPAIRED_READS_EXAMINED',
                  'UNMAPPED_READS',  
                  'UNPAIRED_READ_DUPLICATES',
                  'READ_PAIR_OPTICAL_DUPLICATES',    
                  'ESTIMATED_LIBRARY_SIZE']

    #filter out lines that start with # or are just new lines    
    datalines=filter(lambda x: '#' not in x, lines)
    datalines=filter(lambda x: x!='\n' , datalines)

    #then strip newlines from the lines that contain data
    datalines=map(lambda x: x.strip() , datalines)

    keys=datalines[0].split('\t')
    values=datalines[1].split('\t')

    #metrics_dict=dict(zip(keys,values))
    metrics_dict = dict(izip_longest(keys,values, fillvalue='NA'))

    #Only print the keys we care about:
    output_dict = dict(zip(variant_keys,itemgetter(*variant_keys)(metrics_dict)))
    return output_dict, variant_keys

