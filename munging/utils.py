from collections import Iterable
from datetime import datetime
import re
import os
import shutil
import logging
from collections import namedtuple
from munging.annotation import multi_split

from __init__ import __version__

log = logging.getLogger(__name__)

ASSAYS={'OPX':'OncoPlex',    
        'BRO':'ColoSeq',
        'EPI':'EpiPlex',
        'MRW':'MarrowSeq',
        'IMD':'ImmunoPlex',
        'HP':'HotSpot-Heme',
        'STH':'HotSpot-Heme',
        'GLT':'HotSpot-Hereditary',
        'TESTDATA':'testdata',
        'MSI-PLUS':'msi-plus'}

MACHINE_CODES={'H':'hiseq',
               'M':'miseq',
               'N':'nextseq'}

ASSAY_CODES={'colo':'coloseq',
             'onco':'oncoplex',
             'epi':'epiplex',
             'imm':'immunoplex',
             'marrow':'marrowseq',
             'broca-hr':'broca-hr',
             'hp':'hotspot-heme',
             'glt':'hotspot-hereditary',
             'sth':'hotspot-heme'}


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


def munge_pfx(pfx):
    """
    Get the plate,well, library-version, assay, control 
    and machine-run from the pfx
    """
    output=pfx.split('.')[0].split('_')
    if len(output)==5:
        #New control samples hit this
        keys=['sample_id','well','library-version','control','machine-run']
        pfx_info = dict(zip(keys,output))
        pfx_info['mini-pfx']='{sample_id}_{control}'.format(**pfx_info)
        pfx_info['run']=pfx_info['sample_id'][:-2]
        pfx_info['pfx']='{sample_id}_{well}_{library-version}_{control}_{machine-run}'.format(**pfx_info)
        pfx_info['assay']=[value for key,value in ASSAYS.items() if re.search(key, pfx_info['library-version'])][0]
    elif len(output)==4:
        #New non-control samples hit this
        keys=['sample_id','well','library-version','machine-run']
        pfx_info = dict(zip(keys,output))
        pfx_info['mini-pfx']='{sample_id}'.format(**pfx_info)
        pfx_info['run']=pfx_info['sample_id'][:-2]
        pfx_info['pfx']='{sample_id}_{well}_{library-version}_{machine-run}'.format(**pfx_info)
        pfx_info['assay']=[value for key,value in ASSAYS.items() if re.search(key, pfx_info['library-version'])][0]
    elif len(output)==3:
        #Research non-control samples often hit this
        keys=['sample_id','well','library-version']
        pfx_info = dict(zip(keys,output))
        pfx_info['mini-pfx']='{sample_id}'.format(**pfx_info)
        pfx_info['run']=pfx_info['sample_id'][:-2]
        pfx_info['pfx']='{sample_id}_{well}_{library-version}'.format(**pfx_info)
        pfx_info['assay']=[value for key,value in ASSAYS.items() if re.search(key, pfx_info['library-version'])][0]
    elif len(output)==2:
        #MSI-Plus will hit this
        keys=['sample_id', 'library-version']
        pfx_info = dict(zip(keys,output))
        pfx_info['mini-pfx']='{sample_id}'.format(**pfx_info)
        pfx_info['pfx']='{sample_id}'.format(**pfx_info)
        if re.search('msi',pfx_info['library-version'], re.IGNORECASE):
            pfx_info['assay']=[value for key,value in ASSAYS.items() if re.search(key, pfx_info['library-version'])]
    elif len(output)==1:
        #Only the old LMG/OPX should hit this. 
        pfx=output[0]
        pfx_info={'mini-pfx':pfx,
                  'pfx':pfx,
                  'sample_id':pfx}
        if re.search('LMG', pfx):
            pfx_info['assay']='Coloseq'
        elif re.search('OPX', pfx):
            pfx_info['assay']='OncoPlex'
    else:
        print "pfx:", pfx
        raise ValueError('Incorrect pfx given. Expected Plate_Well_Assay_<CONTROL>_MachinePlate.file-type.file-ext')


    return pfx_info

def munge_date(date):
    """
    Convert new date YYMMDD to YYYY-MM-DD
    """
    #only convert if new date format
    if len(date)<7:
        d = datetime.strptime(date, '%y%m%d')
        return d.strftime('%Y-%m-%d')
    else:
        return date


def munge_old_path(pth):
    """
    Get date, run, project, machine, assay, from old-path
    """
    
def munge_path(pth):
    """
    Get date, run, project, machine, assay, prep-type from path
    """
    output=multi_split(pth, '/_')
    #Assuming we want YYMMDD_RUN_PROJECT
    if output[-1]=='output':
        output=output[-4:-1]
        keys=['date','run', 'project']
    #If old version of data that isn't in a 'output' subfolder
    elif len(output)==5:
        keys=['date','machine','assay','run','version','project']
        #check that the third item is the assay
        if not output[2].lower() in ASSAY_CODES.values():
            #is the second item is the assay?
            if not output[1].lower() in ASSAY_CODES.values():
                print "not a good path", output
            project=output[1]+output[3].strip('run')
            output=[output[0],output[2],output[1],output[3], output[4],project] 
        else:
            project=output[2]+output[3].strip('run')
            output=[output[0],output[1],output[2],output[3], output[4],project] 
    else:
        output=output[-3:]
        keys=['date','run', 'project']
    pathinfo = dict(zip(keys,output))
    pathinfo['date']=munge_date(pathinfo['date'])
    pathinfo['project']=pathinfo['project'].lower()
    #Set Machine
    if re.search('[%s][A-Z]\d+' % ''.join(MACHINE_CODES.keys()), pathinfo['run']):
        pathinfo['machine']= MACHINE_CODES[pathinfo['run'][0]]
    #Set assay
    for a in ASSAY_CODES.keys():
        if re.search(a, pathinfo['project'], re.IGNORECASE):
            pathinfo['assay'] = ASSAY_CODES[a]
        elif re.search('MSI', pathinfo['project'], re.IGNORECASE):
            pathinfo['assay'] = 'msi-plus'
    #Set prep type
    if re.search('kapa', pathinfo['project']):
        pathinfo['prep_type']='kapa'
    elif re.search('hotspot', pathinfo['assay']):
        pathinfo['prep_type']='truseq'
    else:
        pathinfo['prep_type']='sure_select'

    return pathinfo


#############################
def cast(val):
    """Attempt to coerce `val` into a numeric type, or a string stripped
    of whitespace.
    """

    for func in [int, float, lambda x: x.strip(), lambda x: x]:
        try:
            return func(val)
        except ValueError:
            pass

class Opener(object):
    """Factory for creating file objects

    Keyword Arguments:
    - mode -- A string indicating how the file is to be opened. Accepts the
      same values as the builtin open() function.
    - bufsize -- The file's desired buffer size. Accepts the same values as
      the builtin open() function.
    """

    def __init__(self, mode='r', bufsize=-1):
        self._mode = mode
        self._bufsize = bufsize

    def __call__(self, string):
        if string is sys.stdout or string is sys.stdin:
            return string
        elif string == '-':
            return sys.stdin if 'r' in self._mode else sys.stdout
        elif string.endswith('.bz2'):
            return bz2.BZ2File(string, self._mode, self._bufsize)
        elif string.endswith('.gz'):
            return gzip.open(string, self._mode, self._bufsize)
        else:
            return open(string, self._mode, self._bufsize)

    def __repr__(self):
        args = self._mode, self._bufsize
        args_str = ', '.join(repr(arg) for arg in args if arg != -1)
        return '{}({})'.format(type(self).__name__, args_str)


def opener(pth, mode='r', bufsize=-1):
    return Opener(mode, bufsize)(pth)

def make_local_copy(outdir, subdir, fname):
    """Copy fname from package data to outdir/subdir (creating dir if
    necessary), and return the path to the copy of fname relative to
    outdir.
    """
    destdir = os.path.join(outdir, subdir)
    mkdir(destdir)
    shutil.copyfile(fname, os.path.join(destdir, fname))
    return os.path.join(subdir, fname)

def create_static_files(file_names, outdir, subdir):
    '''
    Takes a list of files and copies them into a particular directory.
    '''
    file_list = []
    for i in range(len(file_names)):
        file_list.append(make_local_copy(outdir, subdir, file_names[i]))
    return file_list

def include_file(outdir, fname):
    mkdir(outdir)
    dest = os.path.join(outdir, os.path.basename(fname))
    shutil.copyfile(fname, dest)
    return os.path.join(os.path.basename(outdir), os.path.basename(fname))


def timestamp_now():
    """
    Produce a string with date and time information for a report
    """
    return datetime.now().strftime("%A, %B %d, %Y, %I:%M %p")
