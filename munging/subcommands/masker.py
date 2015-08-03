"""
Create masked output files

Usage:

 munge mask_masker $SAVEPATH/$PFX

"""
from shutil import copyfile
import pprint
import logging
import os
import csv
import re
from munging.utils import walker
from itertools import ifilter
import sys
from munging.filters import any_analysis, maskable
log = logging.getLogger(__name__)

MASK_CODES={
    'BROCA':{'Genes':('BRCA1','BRCA2')},
    'EPIV01':{'Genes':('ALDH7A1','ARX','CDKL5','CHD2','FOXG1',
                       'GABRA1','GABRG2','KCNQ2','KCNQ3','KCNT1',
                       'MECP2','MEF2C','PCDH19','PNKP','PNPO',
                       'PTEN','SCN1A','SCN2A','SCN8A','SLC2A1',
                       'SLC9A6','STXBP1','SYNGAP1','TCF4','UBE3A')},
    'MEGPX':{'Genes':('AKT1','AKT3','PTEN','PIK3CA','PIK3R2')},
    'IMDF01':{'Genes':('ADA','AK2','AP3B1','ATM','BLM','BLNK',
                       'BTK','CARD11','CD3D','CD3E','CD3G',' CD8A',
                       'CD27','CD79A','CD79B','CD247','CHD7','CIITA',
                       'CORO1A','CTLA4','DCLRE1C','DOCK8','EBF1',
                       'FOXN1','FOXP3','GATA2','IGHM','IGLL1','IKBKB',
                       'IKBKG,IKZF1','IL2RA','IL2RG','ILR7','ITK',
                       'JAK3','LCK','LIG4','LRRC8A','LYST','MAGT1',
                       'MALT1','MRE11A','NBN','NFKBIA','NHEJ1','ORAI1',
                       'PIK3R1','PNP','PRF1','PRKDC','PTPRC','RAB27A',
                       'RAG1','RAG2','RFX5','RFXANK','RFXAP,RMRP',
                       'SH2D1A','SP110','STAT1','STAT5B','STIM1','STK4',
                       'STX11','STXBP2','TAP1','TAP2','TAPBP','TBX1',
                       'TRAC','TTC7A','UNC13D','XIAP','ZAP70')},
    'IMDS01':{'Genes':('ADA','AK2','ATM','BLM','CD3D','CD3E','CD3G',
                       'CD8A','CD27','CD247','CHD7','CIITA','CORO1A','CTLA4',
                       'DCLRE1C','DOCK8','FOXN1','FOXP3,GATA2','IKBKB',
                       'IKBKG','IL2RA','IL2RG','ILR7','ITK','JAK3','LCK',
                       'LIG4','MAGT1','MALT1','MRE11A','NBN','NFKBIA',
                       'NHEJ1','ORAI1','PNP','PRKDC','PTPRC','RAG1',
                       'RAG2','RFX5','RFXANK','RFXAP','RMRP','SP110',
                       'STAT1','STAT5B','STIM1','STK4','TAP1','TAP2',
                       'TAPBP','TBX1','TRAC','TTC7A','ZAP70')},
    'IMDB01':{'Genes':('BLNK','BTK','CARD11','EBF1','CD79A','CD79B',
                       'GATA2','IGHM','IGLL1','IKZF1','LIG4','LRRC8A',
                       'MALT1','NHEJ1','PIK3R1','PRKDC','RAG1','RAG2','SH2D1A')},
    'IMDH01':{'Genes':('AP3B1','ITK','LYST','MAGT1','PRF1','RAB27A','SH2D1A',
                       'STX11','STXBP2','UNC13D','XIAP')}
}

def build_parser(parser):
    parser.add_argument('path',
                        help='Path to analysis files')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-o','--order_code', 
                       choices=['BROCA','IMDH01','IMDB01','IMDS01',
                                'IMDF01','MEGPX','EPIV01'],
                       help="Order code for genes that were tested")
    group.add_argument('-g','--gene_list', nargs='+',
                       help="List of genes that were tested, must have spece between names")
    #This script filters:
    # 4_SV_Crest
    # 5_SV_Breakdancer
    # 6_SV_Pindel
    # 10_SNP_Indel
    
    #No filtering needed:
    # 1_QC_Metrics
    # 2_QC_by_Gene
    # 3_QC_by_Exons
    # 7_CNV_Gene
    # 8_CNV_Exon
    # 9_Clinically_Flagged
    # 11_MSI

def mask_file_by_gene(data, genes):
    """
    Create dictionary of analysis file info
    Ensure fname['Gene'] exists
    take in dictionary of data,
    return if Gene is in masking genes list

    """
    output=[]
    gene_headers=('Gene','Gene_1','Gene_2', 'Clinically_Flagged')		    
    for d in data:
        for key, value in d.iteritems():   # iter on both keys and values
            if key in gene_headers:
                #some entries have multiple genes, if any are in masking list, include in output
                sample_genes=set(value.split(',')).intersection(genes)
                if sample_genes:
                    if d not in output:
                        output.append(d)
    
    return output

def action(args):
    infiles = walker(args.path)  
    files=[i for i in infiles]
    if len(files)<1:
        print "No files where found. Are there subfolders for each sample?"
        sys.exit(1)
    #Get the set of genes for masking, based on cli entry
    if args.order_code:
        mask=MASK_CODES[args.order_code]['Genes'] 
    elif args.gene_list:
        mask=args.gene_list 
    print 'Genes in output: %s ' % ([i for i in mask])
    #Grab files for filtering
    files = ifilter(any_analysis, files)

    #Filter to only those that are "maskable"
    files = ifilter(maskable, files)

    for pth in files:
        analysis_type=pth.fname.split('.')[1]
        #In this case, pfx is actually PFX_analysis_type
        pfx = pth.fname.strip('.txt')
        #Create file names for new output
        full_output=os.path.join(pth.dir, (pfx+'.full.txt'))
        masked_output=os.path.join(pth.dir, (pfx+'.masked.txt'))

        #Open data for masking
        data=csv.DictReader(open(os.path.join(pth.dir,pth.fname)),delimiter='\t')

        #Open output for writing
        writer = csv.DictWriter(open(masked_output, 'w'),
                                fieldnames=data.fieldnames,
                                quoting=csv.QUOTE_MINIMAL,
                                extrasaction='ignore',
                                delimiter='\t')
        writer.writeheader()
        # #Mask data
        print "filtering %s" % analysis_type
        writer.writerows(mask_file_by_gene(data, mask))
        #Move the files so the masked is Analysis.txt and the full is labeled
        copyfile(os.path.join(pth.dir,pth.fname),full_output)
        os.rename(masked_output, os.path.join(pth.dir,pth.fname))
