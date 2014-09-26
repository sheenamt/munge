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

def only_analysis(pth):
    """
    True only for pfx_Analysis.{txt,csv}
    """

    return pth.fname.split('_')[-2] not in set(['CNV', 'Genotype', 'Pindel', 'SV', 'QC', 'Quality', 'Breakdancer', 'quality', 'genotype','Gene', 'Exon', 'MSI', 'Ethnicity'])

def cnv_gene_analysis(pth):
    """
    True only for pfx_CNV_Gene_Analysis.{txt,csv}
    """
    return pth.fname.split('_')[-2] in set (['Gene']) and pth.fname.split('_') [-3] not in set(['QC'])

def cnv_exon_analysis(pth):
    """
    True only for pfx_CNV_Exon_Analysis.{txt,csv}
    """
    return pth.fname.split('_')[-2] in set (['Exon']) and pth.fname.split('_') [-3] not in set(['QC'])

def cnv_bins_analysis(pth):
    """
    True only for pfx_CNV_Bins.txt
    """
    return pth.fname.split('_')[-2] in set (['Bins'])

def pindel_analysis(pth):
    """
    True only for pfx_Pindel_Analysis.{txt,csv}
    """
    return pth.fname.split('_')[-2] in set (['Pindel'])

def msi_analysis(pth):
    """
    True only for pfx_MSI_Analysis.{txt,csv}
    """
    return pth.fname.split('_')[-2] in set (['MSI'])


