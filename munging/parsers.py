"""
Parsers for the top level summary files
pindel, cnv_exon, cnv_gene, snp, msi, quality, clin_flagged
"""
import os
import csv
import re
import sys
import copy

from itertools import count, groupby, chain, ifilter , izip_longest
from operator import itemgetter

from munging import filters
from munging.utils import munge_pfx, walker

"""Each function parses a group of sample files for desired information,
grouping based on the variant_keys list,
some include additional annotation headers,
sample counts, and scores calculated based on counts
"""

def parse_quality(files, specimens, annotation, prefixes, variant_keys, sort_order):
    """ Parse the sample quality analysis file, from hs_metrics"""
    files = filter(filters.quality_analysis, files)
    variant_keys = ['MEAN_TARGET_COVERAGE']
    for sample in sort_order:
        #Grab the file for each sample, in specified sort order
        pfx_file = [s for s in files if sample in s.fname]
        if pfx_file:
            pfx_file = pfx_file[0]
            pfx = munge_pfx(pfx_file.fname)
            log_pfx=pfx['mini-pfx']
            prefixes.append(log_pfx)
            with open(os.path.join(pfx_file.dir, pfx_file.fname)) as fname:
                reader = csv.DictReader(fname, delimiter='\t')
                for row in reader:
                    if re.search('version',row['MEAN_TARGET_COVERAGE']):
                        continue
                    variant = tuple(k for k in variant_keys)
                    specimens[variant][log_pfx] = row['MEAN_TARGET_COVERAGE']
                    annotation[variant] = specimens[variant]

    fieldnames = variant_keys + prefixes
    return specimens, annotation, prefixes, fieldnames, variant_keys 

def parse_clin_flagged(files, specimens, annotation, prefixes, variant_keys, sort_order):
    """Parse the Genotype output, which is the reads of clin_flagged found"""
    files = filter(filters.genotype_analysis, files)
    variant_keys = ['Position','Ref_Base','Var_Base' ]
    for sample in sort_order:
        #Grab the file for each sample, in specified sort order
        pfx_file = [s for s in files if sample in s.fname]
        if pfx_file:
            pfx_file = pfx_file[0]
            pfx = munge_pfx(pfx_file.fname)
            #Create a smaller version of this really long string
            reads_pfx=pfx['mini-pfx']+'_Variants'
            prefixes.append(reads_pfx)
            with open(os.path.join(pfx_file.dir, pfx_file.fname)) as fname:
                reader = csv.DictReader(fname, delimiter='\t')
                for row in reader:
                    variant = tuple(row[k] for k in variant_keys)
                    specimens[variant][reads_pfx]=row['Variant_Reads']
                    annotation[variant] = row

    annotation_headers = ['Clinically_Flagged']
    fieldnames = variant_keys + annotation_headers + prefixes
    return specimens, annotation, prefixes, fieldnames, variant_keys            

def parse_msi_flagged(files, specimens, annotation, prefixes, variant_keys, sort_order):
    """Parse the Genotype output, which is the reads of clin_flagged found"""
    files = filter(filters.genotype_analysis, files)
    variant_keys = ['Position','Ref_Base','Var_Base' ]
    for sample in sort_order:
        #Grab the file for each sample, in specified sort order
        pfx_file = [s for s in files if sample in s.fname]
        if pfx_file:
            pfx_file = pfx_file[0]
            pfx = munge_pfx(pfx_file.fname)
            #Create a smaller version of this really long string
            reads_pfx=pfx['mini-pfx']+'_Variants|Total'
            status_pfx=pfx['mini-pfx']+'_Status'
            prefixes.append(reads_pfx)
            prefixes.append(status_pfx)
            with open(os.path.join(pfx_file.dir, pfx_file.fname)) as fname:
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

def parse_hotspot_flagged(files, specimens, annotation, prefixes, variant_keys, sort_order):
    """Parse the Genotype output, which is the reads of clin_flagged found"""
    files = filter(filters.genotype_analysis, files)
    variant_keys = ['Position','Ref_Base','Var_Base' ]
    for sample in sort_order:
        #Grab the file for each sample, in specified sort order
        pfx_file = [s for s in files if sample in s.fname]
        if pfx_file:
            pfx_file = pfx_file[0]
            pfx = munge_pfx(pfx_file.fname)
            #Create a smaller version of this really long string
            reads_pfx=pfx['mini-pfx']+'_Variants|Total'
            status_pfx=pfx['mini-pfx']+'_Status'
            prefixes.append(reads_pfx)
            prefixes.append(status_pfx)
            with open(os.path.join(pfx_file.dir, pfx_file.fname)) as fname:
                reader = csv.DictReader(fname, delimiter='\t')
                for row in reader:
                    variant = tuple(row[k] for k in variant_keys)
                    #Skip lines with no read info
                    if not row['Variant_Reads'] and not row['Valid_Reads']:
                        continue
                    try:
                        frac = "{0:.4f}".format(float(row['Variant_Reads'])/float(row['Valid_Reads']))
                    except ZeroDivisionError:
                        frac = '0'
                    specimens[variant][reads_pfx]=row['Variant_Reads']+'|'+row['Valid_Reads']
                    if int(row['Valid_Reads']) >= 100:
                        if float(frac) > 0.70:
                            specimens[variant][status_pfx]='HOMO'
                        elif float(frac) <= 0.10:
                            specimens[variant][status_pfx]='NEG'
                        elif 0.20 <= float(frac) <= 0.70 :
                            specimens[variant][status_pfx]='HET'
                        else:
                            specimens[variant][status_pfx]='REVIEW'                    
                    else:
                        specimens[variant][status_pfx]='REVIEW'
                    annotation[variant] = row
    annotation_headers = ['Clinically_Flagged']
    fieldnames = variant_keys + annotation_headers + prefixes
    return specimens, annotation, prefixes, fieldnames, variant_keys            

def parse_glt_flagged(files, specimens, annotation, prefixes, variant_keys, sort_order):
    """Parse the Genotype output, which is the reads of clin_flagged found"""

    files = filter(filters.genotype_analysis, files)
    variant_keys = ['Position','Ref_Base','Var_Base' ]
    for sample in sort_order:
        #Grab the file for each sample, in specified sort order
        pfx_file = [s for s in files if sample in s.fname]
        if pfx_file:
            pfx_file = pfx_file[0]
            pfx = munge_pfx(pfx_file.fname)
            #Create a smaller version of this really long string
            reads_pfx=pfx['mini-pfx']+'_Variants|Total'
            status_pfx=pfx['mini-pfx']+'_Status'
            prefixes.append(reads_pfx)
            prefixes.append(status_pfx)
            with open(os.path.join(pfx_file.dir, pfx_file.fname)) as fname:
                reader = csv.DictReader(fname, delimiter='\t')
                for row in reader:
                    variant = tuple(row[k] for k in variant_keys)
                    #Skip lines with no read info
                    if not row['Variant_Reads'] and not row['Valid_Reads']:
                        continue
                    try:
                        frac = "{0:.4f}".format(float(row['Variant_Reads'])/float(row['Valid_Reads']))
                    except ZeroDivisionError:
                        frac = '0'
                    specimens[variant][reads_pfx]=row['Variant_Reads']+'|'+row['Valid_Reads']
                    if int(row['Valid_Reads']) >= 100:
                        if float(frac) >= 0.98:
                            specimens[variant][status_pfx]='HOMO'
                        elif float(frac) <= 0.10:
                            specimens[variant][status_pfx]='NEG'
                        elif 0.40 <= float(frac) <= 0.65 :
                            specimens[variant][status_pfx]='HET'
                        else:
                            specimens[variant][status_pfx]='REVIEW'                    
                    else:
                        specimens[variant][status_pfx]='REVIEW'
                    annotation[variant] = row
    annotation_headers = ['Clinically_Flagged']
    fieldnames = variant_keys + annotation_headers + prefixes
    return specimens, annotation, prefixes, fieldnames, variant_keys            

def parse_pindel(files, specimens, annotation, prefixes, variant_keys, sort_order):
    """Parse the pindel analysis file, give total counts of samples with site"""

    files = filter(filters.pindel_analysis, files)
    variant_keys = ['Position', 'Gene', 'Size']
    #Other annotation to keep 
    annotation_headers = [
        'Gene_Region',
        'Event_Type',
        'Transcripts'
        ]
    for sample in sort_order:
        #Grab the file for each sample, in specified sort order
        pfx_file = [s for s in files if sample in s.fname]
        if pfx_file:
            pfx_file = pfx_file[0]
            pfx = munge_pfx(pfx_file.fname)
            #Create a smaller version of this really long string
            prefixes.append(pfx['mini-pfx'])
            with open(os.path.join(pfx_file.dir, pfx_file.fname)) as fname:
                reader = csv.DictReader(fname, delimiter='\t')
                for row in reader:
                    variant = tuple(row[k] for k in variant_keys)
                    #Update the specimen dict for this variant, for this pfx, report the Reads found
                    try:
                        specimens[variant][pfx['mini-pfx']] = max(row['bbmergedReads'], row['bwamemReads'])
                    except KeyError:
                        specimens[variant][pfx['mini-pfx']] = row['Reads']
                    annotation[variant] = row

    #Update the specimen dict for this variant, count samples present
    for key, value in specimens.iteritems():
        specimens[key]['Count']=len(value)

    #Add 'Count' to prefixes for correct dict zipping/printing    
    prefixes.append('Count')
    fieldnames = variant_keys + annotation_headers + prefixes
    return specimens, annotation, prefixes, fieldnames, variant_keys            

def parse_snp(files, specimens, annotation, prefixes, variant_keys, sort_order):
    """Parse the snp output file, give ref|var read counts per sample"""
    files = filter(filters.snp_analysis, files)
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
        'UW_Freq',
        'UW_Count',
        'UW_DEC_p',
        '1000g_ALL',
        'EVS_esp6500_ALL',
        '1000g_AMR',
        'EVS_esp6500_AA',
        '1000g_EUR',
        'EVS_esp6500_EU',
        '1000g_SAS',
        '1000g_EAS',
        '1000g_AFR',
        'ADA_Alter_Splice',
        'RF_Alter_Splice',
]

    for sample in sort_order:
        #Grab the file for each sample, in specified sort order
        pfx_file = [s for s in files if sample in s.fname]
        if pfx_file:
            pfx_file = pfx_file[0]
            pfx = munge_pfx(pfx_file.fname)
            reads_pfx=pfx['mini-pfx']+'_Ref|Var'
            prefixes.append(reads_pfx)
            with open(os.path.join(pfx_file.dir, pfx_file.fname)) as fname:
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


def parse_cnv_exon(files, specimens, annotation, prefixes, variant_keys, sort_order):
    """Parse the cnv_exon output, give ave_log_ratio"""
    files = filter(filters.cnv_exon_analysis, files)
    variant_keys = ['Position', 'Gene' ]
    for sample in sort_order:
        #Grab the file for each sample, in specified sort order
        pfx_file = [s for s in files if sample in s.fname]
        if pfx_file:
            pfx_file = pfx_file[0]
            pfx = munge_pfx(pfx_file.fname)
            #Create a smaller version of this really long string
            log_pfx=pfx['mini-pfx']+'_Log'
            prefixes.append(log_pfx)
            with open(os.path.join(pfx_file.dir, pfx_file.fname)) as fname:
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

def parse_cnv_gene(files, specimens, annotation, prefixes, variant_keys, sort_order):
    """Parse the cnv_genes output, give ave_log_ratio"""
    files = filter(filters.cnv_gene_analysis, files)
    variant_keys = ['Position', 'Gene' ]
    for sample in sort_order:
        #Grab the file for each sample, in specified sort order
        pfx_file = [s for s in files if sample in s.fname]
        if pfx_file:
            pfx_file = pfx_file[0]
            pfx = munge_pfx(pfx_file.fname)
            #Create a smaller version of this really long string
            log_pfx=pfx['mini-pfx']+'_Log'
            prefixes.append(log_pfx)
            with open(os.path.join(pfx_file.dir, pfx_file.fname)) as fname:
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

def parse_annotsv(files, specimens, annotation, prefixes, variant_keys, sort_order):
    """Parse the annotsv analysis file, give total counts of samples with site"""

    files = filter(filters.annotsv_analysis, files)
    variant_keys = ['Event1', 'Event2', 'Gene1', 'Gene2', 'NM']
    for sample in sort_order:
        #Grab the file for each sample, in specified sort order
        pfx_file = [s for s in files if sample in s.fname]
        if pfx_file:
            pfx_file = pfx_file[0]
            pfx = munge_pfx(pfx_file.fname)
            #Create a smaller version of this really long string
            prefixes.append(pfx['mini-pfx'])
            with open(os.path.join(pfx_file.dir, pfx_file.fname)) as fname:
                reader = csv.DictReader(fname, delimiter='\t')
                for row in reader:
                    variant = tuple(row[k] for k in variant_keys)
                    specimens[variant][pfx['mini-pfx']] = row['QUAL']
                    annotation[variant] = row

                    
    #Update the specimen dict for this variant, count samples present
    for key, value in specimens.iteritems():
        specimens[key]['Count']=len(value)

    #Add 'Count' to prefixes for correct dict zipping/printing    
    prefixes.append('Count')
    fieldnames = variant_keys + prefixes
    return specimens, annotation, prefixes, fieldnames, variant_keys            
