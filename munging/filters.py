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
    True only for pfx_Analysis.{txt,csv}
    """

    return pth.fname.split('.')[-2] in set(['SNP_Analysis'])

def cnv_gene_analysis(pth):
    """
    True only for pfx_CNV_Gene_Analysis.{txt,csv}
    """
    return pth.fname.split('.')[-2] in set (['CNV_Gene_Analysis'])

def cnv_exon_analysis(pth):
    """
    True only for pfx_CNV_Exon_Analysis.{txt,csv}
    """

    return pth.fname.split('.')[-2] in set (['CNV_Exon_Analysis'])

def cnv_bins_analysis(pth):
    """
    True only for pfx_CNV_Bins.txt
    """
    return pth.fname.split('.')[-2] in set (['Bins'])

def pindel_analysis(pth):
    """
    True only for pfx_Pindel_Analysis.{txt,csv}
    """
    return pth.fname.split('.')[-2] in set (['Pindel_Analysis'])

def msi_analysis(pth):
    """
    True only for pfx_MSI_Analysis.{txt,csv}
    """
    return pth.fname.split('.')[-2] in set (['MSI_Analysis'])

def msi_file_finder(pth):
    """
    Return True if pth represents an analysis file.
    """
    return bool(re.search(r'.msi.txt', pth.fname))

def quality_analysis(pth):
    """
    True only for pfx_Quality_Analysis.{txt,csv}
    """

    return pth.fname.split('.')[-2] in set (['Quality_Analysis'])

def genotype_analysis(pth):
    """
    True only for pfx_Genotype_Analysis.{txt,csv}
    """
    return pth.fname.split('.')[-2] in set (['Genotype_Analysis'])


