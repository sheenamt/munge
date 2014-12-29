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

def build_parser(parser):
    parser.add_argument('infiles', action='append', nargs='+',
                        help='Input files')
    
    parser.add_argument('order_code', choices=['BROCA'],
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
    
    #No filtering:
    # 1_QC_Metrics
    # 11_MSI

def mask_file_by_gene(fname, genes, outfile):
    """
    Create dictionary of analysis file info
    Ensure fname['Gene'] exists
    take in dictionary of data,
    return if Gene is in masking genes list

    """
    data=csv.DictReader(open(fname), delimiter='\t')
    gene_headers=('Gene','Gene_1','Gene_2', 'Clinically_Flagged')
    print 'outfile:', outfile
    out_data=[]
    output = csv.DictWriter(open(outfile, 'w'),
                            fieldnames=data.fieldnames,
                            quoting=csv.QUOTE_MINIMAL,
                            extrasaction='ignore',
                            delimiter='\t')
    output.writeheader()
    for d in data:
        for key, value in d.iteritems():   # iter on both keys and values
            if key in gene_headers and value in genes:
                out_data.append(d)
                for d in out_data:
                    output.writerow(d)

    def action(args):
        (infiles, ) = args.infiles
        ordered_genes={'BROCA':{'Genes':('BRCA1','BRCA2')}}
        mask=ordered_genes[args.order_code]['Genes']
        print 'Genes in output: %s ' % ([i for i in mask])
        
        for fname in infiles:
            (f_path, f_name) = os.path.split(fname)
            if re.search('Analysis', f_name):
                print 'parsing:', f_name
                (analysis_type,ext) = os.path.splitext(f_name)
                #Rename Analysis file to Analysis_full
                full_input=os.path.join(f_path, (analysis_type+'_full'+ext))
                masked_output=os.path.join(f_path, (analysis_type+'_masked'+ext))
                #Mask data
                mask_file_by_gene(fname, mask, masked_output)
                #Move the files so the masked is Analysis.txt and the full is labeled
                copyfile(fname,full_input)
                os.rename(masked_output, fname)
