"""
Create xlsx workbook from all output files

usage:

 munge xlsmaker (Analysis/Combined) $SAVEPATH/$PFX* -o $SAVEPATH/${PFX}_Analysis.xls

"""
import csv
import re
import os

from xlsxwriter import Workbook
from munging.annotation import build_variant_id, multi_split

def build_parser(parser):
    parser.add_argument(
        'type', choices=['Analysis', 'Combined'],
        help='Type of files to create workbook, Analysis or Combined')
    parser.add_argument(
        'infiles', action='append', nargs='+',
        help='Input files')
    parser.add_argument(
        '-o', '--outfile',
        help='Output file', required=True)


def float_if_possible(strg):
    """
    Convert integers to float instead of string where applicable.
    """
    strg=strg.decode('utf-8')
    try:
        value=float(strg)
    except ValueError:
        value=strg
    return value

def process_files(infiles, tab, filetype):
    """
    Rename the analysis files for workbook
    """
    for fname in infiles:
        (f_path, f_name) = os.path.split(fname)
        if re.search(str(filetype), f_name):
            (f_short_name, f_extension) = os.path.splitext(f_name)
            sheet_name = multi_split(f_short_name, '._')
            try:
                #48_A03_BROv7_HA0186_NA12878_QC_Analysis
                if sheet_name[-2] == 'QC':
                    sheet_name = '0_QC'
                #48_A03_BROv7_HA0186_NA12878_Quality_Analysis
                elif sheet_name[-2] == 'Quality':
                    sheet_name = '1_QC_Metrics'
                #Gene/Exon is either coverage kit or CNV, can't be both
                #NA12878-HP998-HHv2.Coverage_Gene
                elif sheet_name[-3]=='Coverage':
                    if sheet_name[-2]=='Gene':
                        sheet_name = '2_QC_by_Gene'
                    elif sheet_name[-2]=='Exon':
                        sheet_name = '3_QC_by_Exon'
                #48_A03_BROv7_HA0186_NA12878_CNV_QC_[Exon/Gene]_Analysis
                elif sheet_name[-3] == 'QC':
                    if sheet_name[-2] == 'Gene':
                        sheet_name = '2_QC_by_Gene'
                    elif sheet_name[-2] == 'Exon':
                        sheet_name = '3_QC_by_Exon'
                #48_A03_BROv7_HA0186_NA12878_CNV_[Exon/Gene]_Analysis
                elif sheet_name[-3] == 'CNV':
                    if sheet_name[-2] == 'Gene':
                        sheet_name = '7_CNV_Gene'
                    elif sheet_name[-2] == 'Exon':
                        sheet_name = '8_CNV_Exon'
                #48_A03_BROv7_HA0186_NA12878_SV_Analysis
                elif sheet_name[-2] == 'SV':
                    sheet_name = '4_SV_Gridss'
                #48_A03_BROv7_HA0186_NA12878_Breakdancer_Analysis
                elif sheet_name[-2] == 'Breakdancer':
                    sheet_name = '5_SV_Breakdancer'
                #48_A03_BROv7_HA0186_NA12878_Pindel_Analysis
                elif sheet_name[-2] == 'Pindel':
                    sheet_name = '6_SV_Pindel'
                #48_A03_BROv7_HA0186_NA12878_Genotype_Analysis
                elif sheet_name[-2] == 'Genotype':
                    sheet_name = '9_Clinically_Flagged'
                #48_A03_BROv7_HA0186_NA12878_MSI_Analysis
                elif sheet_name[-2] == 'MSI':
                    sheet_name = '12_MSI'
                #48_A03_BROv7_HA0186_NA12878_Amplicon_Analysis
                elif sheet_name[-2] == 'Amplicon':
                    sheet_name = '13_Amplicons'
                #48_A03_BROv7_HA0186_NA12878_PolyHunter_Analysis
                elif sheet_name[-2] == 'PolyHunter':
                    sheet_name = '14_PolyHunter'
                #48_A03_BROv7_HA0186_NA12878_SNP_Indel_Analysis
                elif sheet_name[-2] == 'SNP':
                    sheet_name = '10_SNP'
                #48_A03_BROv7_HA0186_NA12878_INDEL_Analysis
                elif sheet_name[-2] == 'INDEL':
                    sheet_name = '11_INDEL'
                if sheet_name == tab:
                    return sheet_name, fname

            except IndexError:
                continue

def add_links(Reader, sheet, fname):
    """
    Process analysis file to add link columns
    """
    #Need to get the pfx from the fname
    (f_path, f_name) = os.path.split(fname)
    f_short_name = f_name.replace('.SNP_Analysis.txt','')

    #rowx is the row#, row is the data
    for rowx, row in enumerate(Reader):
        variant_id, ref_reads, var_reads=build_variant_id(row)
        row.insert(0, variant_id)
        for colx, value in enumerate(row):
            if colx == 0 and not value == 'gendb_link':
                if len(value) < 198:
                    sheet.write_formula(rowx, colx, 'HYPERLINK("https://gendb.labmed.uw.edu/genetics_db/search?variant_id={}","gendb_link")'.format(value))
            elif colx == 15 and not value == 'Faves_Y/N':
                if len(value) < 198:
                    sheet.write_formula(rowx, colx, 'HYPERLINK("https://control.labmed.uw.edu/var_clin/fav/?variant_id={}&pfx={}&ref_reads={}&var_reads={}","mark_fav")'.format(variant_id, f_short_name, ref_reads, var_reads))

            else:
                sheet.write(rowx, colx, float_if_possible(value))

def variant_id_link(Reader, sheet):
    """
    Process analysis file to add link column
    """    
    for rowx, row in enumerate(Reader):
        row.insert(0, build_variant_id(row))
        for colx, value in enumerate(row):
            if colx == 0 and not value == 'link':
                if len(value) < 198:
                    sheet.write_formula(rowx, colx, 'HYPERLINK("https://apps.labmed.uw.edu/genetics_db/search?variant_id={}","link")'.format(value))
            else:
                sheet.write(rowx, colx, float_if_possible(value))

def write_workbook(sheet_name, book, fname):
    """
    Write analysis file as sheet in workbook
    """
    sheet = book.add_worksheet(sheet_name)    
    Reader = csv.reader(open(fname, 'rU'), delimiter='\t')
    if sheet_name == '10_SNP':
        Reader = add_links(Reader, sheet, fname)
    else:
        for rowx, row in enumerate(Reader):
            for colx, value in enumerate(row):
                sheet.write(rowx, colx, float_if_possible(value))

def action(args):

    book = Workbook()
    filetype = args.type
    (infiles, ) = args.infiles
    if filetype == 'Analysis':
        tabs = ['0_QC', '1_QC_Metrics', '2_QC_by_Gene', '3_QC_by_Exon',
                '4_SV_Gridss', '5_SV_Breakdancer', '6_SV_Pindel',
                '7_CNV_Gene', '8_CNV_Exon', '9_Clinically_Flagged', 
                '10_SNP','11_INDEL', '12_MSI', '13_Amplicons', '14_PolyHunter']
        #    for each tab, find its file and process.
        for tab in tabs:
            try:
                #Find file in infiles
                sheet_name, fname = process_files(infiles, tab, filetype)
                print sheet_name, fname
                write_workbook(sheet_name, book, fname)
            except TypeError:
                print "Tab %s not processed" % tab
    elif filetype == 'Combined':
        for fname in infiles:
            (f_path, f_name) = os.path.split(fname)
            if re.search(str(filetype), f_name):
                (f_short_name, f_extension) = os.path.splitext(f_name)
                sheet_name = f_short_name.split('Combined_')
                sheet_name = '_'.join(sheet_name[1:30])
                print sheet_name, fname
                write_workbook(sheet_name, book, fname)
    book.filename=args.outfile

    book.close()
