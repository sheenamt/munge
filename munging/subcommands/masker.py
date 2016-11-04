"""
Create masked output files

Usage:

munge mask_masker file-to-mask order-code/gene-list

"""
import argparse
import sys
import logging
import os
import csv
from munging.utils import walker, validate_gene_list
from itertools import ifilter
from shutil import copyfile
from munging.filters import any_analysis, maskable

log = logging.getLogger(__name__)

MASK_CODES={
    'EPIPX1':('ALDH5A1', 'ALDH7A1', 'ALG13', 'ARHGEF9', 'ARX', 
              'CACNA1A', 'CDKL5', 'CHD2', 'CHRNA2', 'CHRNA4', 
              'CHRNA7', 'CHRNB2', 'CYFIP1', 'DEPDC5', 'DNM1', 
              'EEF1A2', 'FOXG1', 'GABRA1', 'GABRB1', 'GABRB3',
              'GABRG2', 'GNAO1', 'GRIN1', 'GRIN2A', 'GRIN2B', 
              'HCN1', 'HNRNPU', 'IQSEC2', 'KANSL1', 'KCNA2', 
              'KCNB1', 'KCNH1', 'KCNH5', 'KCNQ2', 'KCNQ3', 
              'KCNT1', 'LGI1', 'MBD5', 'MECP2', 'MEF2C', 
              'MTOR', 'NDE1', 'PCDH19', 'PIGA', 'PLCB1', 'PNKP', 
              'PNPO', 'POLG', 'PTEN', 'PURA', 'SCN1A', 'SCN1B', 
              'SCN2A', 'SCN8A', 'SIK1', 'SLC13A5', 'SLC1A2', 
              'SLC25A22', 'SLC2A1', 'SLC35A2', 'SLC6A1', 'SLC9A6', 
              'SPTAN1', 'STX1B', 'STXBP1', 'SYN1', 'SYNGAP1', 
              'TBC1D24', 'TCF4', 'TSC1', 'TSC2', 'UBE3A', 
              'WDR45', 'WWOX', 'ZEB2'),
    'MEGPX1':('ABCC9', 'AKT1', 'AKT2', 'AKT3', 'BRWD3', 
              'CCND2', 'CDKN1C', 'CUL4B', 'DEPDC5', 'DNMT3A', 
              'EED', 'EZH2', 'GLI3', 'GNAQ', 'GNAS', 'GPC3', 
              'HEPACAM', 'KCNJ8', 'KIF7', 'MED12', 'MLC1', 
              'MTOR', 'NFIA', 'NFIX', 'NSD1', 'PIK3CA', 'PIK3R2', 
              'PTCH1', 'PTEN', 'RAB39B', 'RIN2', 'RNF135', 
              'SETD2', 'STRADA', 'TBC1D7', 'TSC1', 'TSC2'),
    'IMDFB1':('ADA','AK2','AP3B1','ATM','BLM','BLNK',
              'BTK','CARD11','CD3D','CD3E','CD3G','CD8A',
              'CD27','CD79A','CD79B','CD247','CHD7','CIITA',
              'CORO1A','CTLA4','DCLRE1C','DOCK8','EBF1',
              'FOXN1','FOXP3','GATA2','IGLL1','IKBKB',
              'IKBKG','IKZF1','IL2RA','IL2RG','IL7R','ITK',
              'JAK3','LCK','LIG4','LRRC8A','LYST','MAGT1',
              'MALT1','MRE11A','NBN','NFKBIA','NHEJ1','ORAI1',
              'PIK3R1','PNP','PRF1','PRKDC','PTPRC','RAB27A',
              'RAG1','RAG2','RFX5','RFXANK','RFXAP','RMRP',
              'SH2D1A','SP110','STAT1','STAT5B','STIM1','STK4',
              'STX11','STXBP2','TAP1','TAP2','TAPBP','TBX1',
              'TTC7A','UNC13D','XIAP','ZAP70'),
    'IMDSB1':('ADA','AK2','ATM','BLM','CD3D','CD3E','CD3G',
              'CD8A','CD27','CD247','CHD7','CIITA','CORO1A','CTLA4',
              'DCLRE1C','DOCK8','FOXN1','FOXP3','GATA2','IKBKB',
              'IKBKG','IL2RA','IL2RG','IL7R','ITK','JAK3','LCK',
              'LIG4','MAGT1','MALT1','MRE11A','NBN','NFKBIA',
              'NHEJ1','ORAI1','PNP','PRKDC','PTPRC','RAG1',
              'RAG2','RFX5','RFXANK','RFXAP','RMRP','SP110',
              'STAT1','STAT5B','STIM1','STK4','TAP1','TAP2',
              'TAPBP','TBX1','TTC7A','ZAP70'),
    'IMDBB1':('BLNK','BTK','CARD11','EBF1','CD79A','CD79B',
              'GATA2','IGLL1','IKZF1','LIG4','LRRC8A',
              'MALT1','NHEJ1','PIK3R1','PRKDC','RAG1','RAG2','SH2D1A'),
    'IMDHB1':('AP3B1','ITK','LYST','MAGT1','PRF1','RAB27A','SH2D1A',
              'STX11','STXBP2','UNC13D','XIAP')
}

def build_parser(parser):
    parser.add_argument(
        '-i', '--infile',
        help='Input file')
    parser.add_argument(
        '-o', '--outfile',
        help='Output file', default=sys.stdout,
        type=argparse.FileType('w'))
    parser.add_argument(
        '-m', '--mask_list', nargs='+',
        help="Order code or list of genes that were tested, must have spece between names")
    parser.add_argument(
        '-g', '--gene_list', 
        help="Tab delimited file of valid genes for this assay, with 'Gene' as header.")
    parser.add_argument(
        '--mask_codes', choices=['1','2'],
        help='Print the built in masking codes (1) and their gene lists (2)')


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
    #If just checking mask code and lists
    if args.mask_codes == '1':
        print 'Codes are: ', MASK_CODES.keys() 
        sys.exit()
    if args.mask_codes == '2':
        for key, value in MASK_CODES.items():
            print 'Code: %s is gene list: \n %s' % (key, value)
        sys.exit()

    infile = args.infile

    if len(infile)<1:
        print "No files where found. Are there subfolders for each sample?"
        sys.exit(1)

    #Get the set of genes for masking, based on cli entry
    try:
        mask=MASK_CODES[args.mask_list[0]]
    except KeyError:
        mask=args.mask_list 

    genes = csv.DictReader(open(args.gene_list,'rU'), delimiter = '\t')
    valid_genes= [i['Gene'] for i in genes]

    #Check that mask code is valide 
    print 'Validating gene list: {}'.format(mask)
    validate_gene_list(mask, valid_genes)

    #Create file names for new output
    full_output=infile
    masked_output=args.outfile
    data=csv.DictReader(open(full_output),delimiter='\t')

    #Open output for writing
    writer = csv.DictWriter(masked_output,
                            fieldnames=data.fieldnames,
                            quoting=csv.QUOTE_MINIMAL,
                            extrasaction='ignore',
                            delimiter='\t')
    writer.writeheader()
    
    # #Mask data
    writer.writerows(mask_file_by_gene(data, mask))
