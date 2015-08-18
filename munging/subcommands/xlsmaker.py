"""
Create xls workbook from all output files

usage:

 munge xlsmaker (Analysis/Combined) $SAVEPATH/$PFX* -o $SAVEPATH/${PFX}_Analysis.xls

"""
import csv
import argparse
import sys
import re
import os

from xlwt import Workbook, Formula
from munging.annotation import build_variant_id, multi_split

book = Workbook()


def build_parser(parser):
    parser.add_argument(
        'type', choices=['Analysis', 'Combined'],
        help='Type of files to create workbook, Analysis or Combined')
    parser.add_argument(
        'infiles', action='append', nargs='+',
        help='Input files')
    parser.add_argument(
        '-o', '--outfile',
        help='Output file', default=sys.stdout,
        type=argparse.FileType('w'))


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
                    sheet_name = '4_SV_Crest'
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
                    sheet_name = '11_MSI'
                #48_A03_BROv7_HA0186_NA12878_Analysis
                elif sheet_name[-2] == 'SNP':
                    sheet_name = '10_SNP_Indel'
                if sheet_name == tab:
                    return sheet_name, fname

            except IndexError:
                continue
def variant_id_link(Reader, sheet):
    """
    Process analysis file to add link column
    """

    for rowx, row in enumerate(Reader):
        row.insert(0, build_variant_id(row))
        for colx, value in enumerate(row):
            if colx == 0 and not value == 'link':
                if len(value) < 198:
                    sheet.write(rowx, colx, Formula('HYPERLINK("https://apps.labmed.uw.edu/genetics_db/search?variant_id={}","link")'.format(value)))
            else:
                sheet.write(rowx, colx, float_if_possible(value))

def write_workbook(sheet_name, fname):
    """
    Write analysis file as sheet in workbook
    """
    sheet = book.add_sheet(sheet_name)
    Reader = csv.reader(open(fname, 'rU'), delimiter='\t')
    if sheet_name == '10_SNP_Indel':
        Reader = variant_id_link(Reader, sheet)
    else:
        for rowx, row in enumerate(Reader):
            for colx, value in enumerate(row):
                sheet.write(rowx, colx, float_if_possible(value))


def action(args):

    filetype = args.type
    (infiles, ) = args.infiles
    if filetype == 'Analysis':
        tabs = ['0_QC', '1_QC_Metrics', '2_QC_by_Gene', '3_QC_by_Exon',
                '4_SV_Crest', '5_SV_Breakdancer', '6_SV_Pindel',
                '7_CNV_Gene', '8_CNV_Exon', '9_Clinically_Flagged', '10_SNP_Indel', '11_MSI']
        #    for each tab, find its file and process.
        for tab in tabs:
            try:
                #Find file in infiles
                sheet_name, fname = process_files(infiles, tab, filetype)
                print sheet_name, fname
                write_workbook(sheet_name, fname)
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
                write_workbook(sheet_name, fname)
    book.save(args.outfile)
