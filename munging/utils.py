from collections import Iterable
import os
import shutil
import logging
from collections import namedtuple
from munging.annotation import multi_split

from __init__ import __version__

log = logging.getLogger(__name__)

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

def munge_path(pth):
    """
    Get run, machine, assay and capture from path
    """
    output=multi_split(pth, '/_')
    if len(output)<4:
        raise ValueError('Incorrect path given. Must be in the format of YYYY-MM-DD_Machine_Assay_Run#_version')
    pathinfo={}
    run=""
    for i in output:
        i=i.lower()
        if i.startswith('20'):
            run=i
        elif i.startswith('run'):
            pathinfo['run']=run+'_'+i
        elif i.startswith('hi') or i.startswith('mi'):
            pathinfo['machine']=i
        elif i.startswith('onco') or i.startswith('colo'):
            pathinfo['assay']=i
        elif i.startswith('v'):
            pathinfo['capture']=i

    return pathinfo

def munge_pfx(pfx):
    """
    Change the pfx output in files to a shorter version
    """
    output=multi_split(pfx, '/_.')
    keys=['sample_id','well','library-version','control','machine-run']
    pfx_info=dict(zip(keys,output))
    pfx_info['control']=check_control(pfx_info['control'])
    pfx_info['mini-pfx']=('_'.join([pfx_info['sample_id'],pfx_info['control']])).strip('_')
    pfx_info['pfx']='_'.join([pfx_info['sample_id'],pfx_info['well'],pfx_info['control'],pfx_info['library-version']])
    pfx_info['run']='{sample_id}'.format(**pfx_info)[:-2]
    return pfx_info

def check_control(control):
    """Check the control is actually the control,
    return it formatted, otherwise return empty string"""

    if control in set(['NA12878']):
        return control
    else:
        return ''
