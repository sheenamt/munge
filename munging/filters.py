"""
Filters for analysis files. In each function, pth is a namedtuple
with attributes (dir, fname).
"""

import re

def any_analysis(pth):
    """
    Return True if pth represents an analysis file.
    """

    return bool(re.search(r'Analysis\.(csv|txt)$', pth.fname))

def snp_analysis(pth):
    """
    True only for pfx.Analysis.{txt,csv}
    """

    return pth.fname.split('.')[-2] in set(['SNP_Analysis'])

def cnv_gene_analysis(pth):
    """
    True only for pfx.CNV_Gene_Analysis.{txt,csv}
    """
    return pth.fname.split('.')[-2] in set (['CNV_Gene_Analysis'])

def cnv_exon_analysis(pth):
    """
    True only for pfx.CNV_Exon_Analysis.{txt,csv}
    """

    return pth.fname.split('.')[-2] in set (['CNV_Exon_Analysis'])

def cnv_bins(pth):
    """
    True only for pfx.CNV_bins.txt
    """
    return pth.fname.split('.')[-2] in set (['CNV_bins'])

def pindel_analysis(pth):
    """
    True only for pfx.Pindel_Analysis.{txt,csv}
    """
    return pth.fname.split('.')[-2] in set (['Pindel_Analysis'])

def msi_analysis(pth):
    """
    True only for pfx.MSI_Analysis.{txt,csv}
    """
    return pth.fname.split('.')[-2] in set (['MSI_Analysis'])

def quality_analysis(pth):
    """
    True only for pfx.Quality_Analysis.{txt,csv}
    """
    return pth.fname.split('.')[-2] in set (['Quality_Analysis'])

def msi_file_finder(pth):
    """
    Return True if pth represents a msi file file.
    """
    return bool(re.search(r'.msi.txt', pth.fname))

def hs_file_finder(pth):
    """
    Return True if pth represents an hs_metrics file.
    """
    return bool(re.search(r'.hs_metrics', pth.fname))

def quality_file_finder(pth):
    """
    Return True if pth represents an hs_metrics file.
    """
    return bool(re.search(r'.quality_metrics', pth.fname))

def genotype_analysis(pth):
    """
    True only for pfx.Genotype_Analysis.{txt,csv}
    """
    return pth.fname.split('.')[-2] in set (['Genotype_Analysis'])


def maskable(pth):
    """
    True only for maskable files
    """
    #This script filters:
    # 4_SV_Crest
    # 5_SV_Breakdancer
    # 6_SV_Pindel
    # 7_CNV_Gene
    # 8_CNV_Exon
    # 10_SNP_Indel

    if pth.fname.split('.')[-2] in set (['SNP_Analysis',
                                         'CNV_Exon_Analysis',
                                         'CNV_Gene_Analysis',
                                         'Pindel_Analysis',
                                         'Breakdancer_Analysis',
                                         'SV_Analysis']):
        return pth.fname.split('.')[-2]


