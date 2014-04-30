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

    return pth.fname.split('_')[-2] not in set(['CNV', 'Genotype', 'Pindel', 'SV', 'QC', 'Quality', 'Breakdancer', 'quality', 'genotype','Gene', 'Exon'])

def cnv_analysis(pth):
    """
    True only for pfx_CNV_Gene_Analysis.{txt,csv}
    """
    return pth.fname.split('_')[-2] in set (['Gene'])

def pindel_analysis(pth):
    """
    True only for pfx_Pindel_Analysis.{txt,csv}
    """
    return pth.fname.split('_')[-2] in set (['Pindel'])


