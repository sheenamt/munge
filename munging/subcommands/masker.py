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


log = logging.getLogger(__name__)

MASK_CODES={
    'BROCA':{'Genes':('BRCA1','BRCA2')},
    'EPIV01':{'Genes':('ALDH7A1','ARX','CDKL5','CHD2','FOXG1',
                       'GABRA1','GABRG2','KCNQ2','KCNQ3','KCNT1',
                       'MECP2','MEF2C','PCDH19','PNKP','PNPO',
                       'PTEN','SCN1A','SCN2A','SCN8A','SLC2A1',
                       'SLC9A6','STXBP1','SYNGAP1','TCF4','UBE3A')},
    'MEGV01':{'Genes':('AKT1','AKT3','PTEN','PIK3CA','PIK3R2')},
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
    parser.add_argument('infiles', action='append', nargs='+',
                        help='Input files')
    
    parser.add_argument('order_code', choices=['BROCA','IMDH01','IMDB01','IMDS01',
                                               'IMDF01','MEGV01','EPIV01'],
                        help="Order code for genes that were tested")

    #This script filters:
    # 4_SV_Crest
    # 5_SV_Breakdancer
    # 6_SV_Pindel
    # 7_CNV_Gene
    # 8_CNV_Exon
    # 10_SNP_Indel
    
    #Not sure for filtering:
    # 9_Clinically_Flagged
    # 2_QC_by_Gene
    # 3_QC_by_Exons
    
    #No filtering needed:
    # 1_QC_Metrics
    # 11_MSI

def mask_file_by_gene(data, genes,out_data):
    """
    Create dictionary of analysis file info
    Ensure fname['Gene'] exists
    take in dictionary of data,
    return if Gene is in masking genes list

    """
    gene_headers=('Gene','Gene_1','Gene_2', 'Clinically_Flagged')		    
    for d in data:
        for key, value in d.iteritems():   # iter on both keys and values
            if key in gene_headers and value in genes:
                out_data.append(d)
    return out_data
                

def action(args):
    (infiles, ) = args.infiles                   
    mask=MASK_CODES[args.order_code]['Genes']
    print 'Genes in output: %s ' % ([i for i in mask])
    out_data=[]

    for fname in infiles:
        (f_path, f_name) = os.path.split(fname)
        if re.search('Analysis', f_name):
            print 'parsing:', f_name
            (analysis_type,ext) = os.path.splitext(f_name)
            #Rename Analysis file to Analysis_full
            full_input=os.path.join(f_path, (analysis_type+'_full'+ext))
            masked_output=os.path.join(f_path, (analysis_type+'_masked'+ext))
            #Mask data
            data=csv.DictReader(open(fname), delimiter='\t')
            out_data=mask_file_by_gene(data, mask,out_data)
            #Write output
            output = csv.DictWriter(open(masked_output, 'w'),
                            fieldnames=data.fieldnames,
                            quoting=csv.QUOTE_MINIMAL,
                            extrasaction='ignore',
                            delimiter='\t')
            output.writeheader()
            output.writerows(out_data)
            
			#Move the files so the masked is Analysis.txt and the full is labeled
            copyfile(fname,full_input)
            os.rename(masked_output, fname)
