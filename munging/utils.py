from collections import Iterable
import datetime
import re
import os
import shutil
import logging
from collections import namedtuple
from munging.annotation import multi_split

from __init__ import __version__

log = logging.getLogger(__name__)

ASSAYS = {'BROv7':'coloseq',
          'BROv8':'coloseq',
          'OPXv3':'oncoplex',
          'OPXv4':'oncoplex',
          'EPIv1':'epiplex',
          'MRWv3':'marrowseq',
          'IMMv1':'immunoplex',
          'MSI-PLUS':'msi-plus'}

def dict_factory(cursor, row):
    """
    Factory to return dicts from sqlite3 queries

    example:

    con = sqlite3.connect(":memory:")
    con.row_factory = dict_factory
    cur = con.cursor()
    cur.execute("select 1 as a")
    print cur.fetchone()["a"]
    """
    return dict((col[0],row[idx]) for idx, col in enumerate(cursor.description))

def flatten(seq):
    """
    Poached from http://stackoverflow.com/questions/2158395/flatten-an-irregular-list-of-lists-in-python

    Don't flatten strings or dict-like objects.
    """
    for el in seq:
        if isinstance(el, Iterable) and not (isinstance(el, basestring) or hasattr(el, 'get')):
            for sub in flatten(el):
                yield sub
        else:
            yield el

def mkdir(dirpath, clobber = False):
    """
    Create a (potentially existing) directory without errors. Raise
    OSError if directory can't be created. If clobber is True, remove
    dirpath if it exists.
    """

    if clobber:
        shutil.rmtree(dirpath, ignore_errors = True)

    try:
        os.mkdir(dirpath)
    except OSError, msg:
        pass

    if not os.path.exists(dirpath):
        raise OSError('Failed to create %s' % dirpath)

    return dirpath

Path = namedtuple('Path', ['dir','fname'])

def walker(dir):
    """Recursively traverse direcory `dir`. For each tuple containing
    (path, dirs, files) return a named tuple with attributes (dir,
    fname) for each file name in `files`.
    """

    for (pth, dirs, files) in os.walk(dir):
        for fname in files:
            yield Path(pth, fname)

def check_control(control):
    """Check the control is actually the control,
    return it formatted, otherwise return empty string"""

    if control in set(['NA12878']):
        return control
    else:
        return None

def munge_pfx(pfx):
    """
    Get the plate,well, library-version, assay, control 
    and machine-run from the pfx
    """
    output=pfx.split('.')[0].split('_')
    if len(output)==5:
        keys=['sample_id','well','library-version','control','machine-run']
        pfx_info = dict(zip(keys,output))
        pfx_info['mini-pfx']='{sample_id}_{control}'.format(**pfx_info)
        pfx_info['run']=pfx_info['sample_id'][:-2]
        pfx_info['pfx']='{sample_id}_{well}_{library-version}_{control}_{machine-run}'.format(**pfx_info)

    elif len(output)==4:
        keys=['sample_id','well','library-version','machine-run']
        pfx_info = dict(zip(keys,output))
        pfx_info['mini-pfx']='{sample_id}'.format(**pfx_info)
        pfx_info['run']=pfx_info['sample_id'][:-2]
        pfx_info['pfx']='{sample_id}_{well}_{library-version}_{machine-run}'.format(**pfx_info)

    else:
        raise ValueError('Incorrect pfx given. Expected Plate_Well_Assay_<CONTROL>_MachinePlate.file-type.file-ext')

    pfx_info['assay']=ASSAYS[pfx_info['library-version']]
    return pfx_info

def munge_date(date):
    """
    Convert new date YYMMDD to YYYY-MM-DD
    """
    #only convert if new date format
    if len(date)<7:
        d = datetime.datetime.strptime(date, '%y%m%d')
        return d.strftime('%Y-%m-%d')
    else:
        return date

def munge_path(pth):
    """
    Get date, run, project, machine, assay, prep-type from path
    """
    output=multi_split(pth, '/_')
    #Assuming we want YYMMDD_RUN_PROJECT
    output=output[-3:]
    keys=['date','run', 'project']
    pathinfo = dict(zip(keys,output))
    pathinfo['date']=munge_date(pathinfo['date'])
    #Set Machine
    if re.search('HA', pathinfo['run']):
        pathinfo['machine']='hiseq'
    elif re.search('MA', pathinfo['run']):
        pathinfo['machine']='miseq'
    #Set assay
    if re.search('colo', pathinfo['project'].lower()):
        pathinfo['assay']='coloseq'
    elif re.search('onco', pathinfo['project'].lower()):
        pathinfo['assay']='oncoplex'
    elif re.search('epi', pathinfo['project'].lower()):
        pathinfo['assay']='epiplex'
    elif re.search('imm', pathinfo['project'].lower()):
        pathinfo['assay']='immunoplex'
    elif re.search('mrw', pathinfo['project'].lower()):
        pathinfo['assay']='marrowseq'
    #Set prep type
    if re.search('kapa', pathinfo['project'].lower()):
        pathinfo['prep_type']='kapa'
    else:
        pathinfo['prep_type']='sure_select'

    return pathinfo
