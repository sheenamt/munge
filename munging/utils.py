from collections import Iterable
import datetime
import re
import os
import shutil
import logging
from collections import namedtuple
import csv
from munging.annotation import multi_split

from __init__ import __version__

log = logging.getLogger(__name__)

ASSAYS = {'BROv7':'coloseq',
          'BROv8':'coloseq',
          'OPXv3':'oncoplex',
          'OPXv4':'oncoplex',
          'OPXv5':'oncoplex',
          'EPIv1':'epiplex',
          'TRU':'truseq',
          'MRWv3':'marrowseq',
          'IMDv1':'immunoplex',
          'TESTDATA':'testdata',
          'MSI-PLUS':'msi-plus'}

MACHINES = {'H':'hiseq',
            'M':'miseq',
            'N':'nextseq'}

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
        pfx_info['assay']=ASSAYS[pfx_info['library-version']]
    elif len(output)==4:
        #New non-control samples hit this
        keys=['sample_id','well','library-version','machine-run']
        pfx_info = dict(zip(keys,output))
        pfx_info['mini-pfx']='{sample_id}'.format(**pfx_info)
        pfx_info['run']=pfx_info['sample_id'][:-2]
        pfx_info['pfx']='{sample_id}_{well}_{library-version}_{machine-run}'.format(**pfx_info)
        pfx_info['assay']=ASSAYS[pfx_info['library-version']]
    elif len(output)==3:
        #Research non-control samples often hit this
        keys=['sample_id','well','library-version']
        pfx_info = dict(zip(keys,output))
        pfx_info['mini-pfx']='{sample_id}'.format(**pfx_info)
        pfx_info['run']=pfx_info['sample_id'][:-2]
        pfx_info['pfx']='{sample_id}_{well}_{library-version}'.format(**pfx_info)
        pfx_info['assay']=ASSAYS[pfx_info['library-version']]
    elif len(output)==2:
        #MSI-Plus will hit this
        keys=['sample_id', 'library-version']
        pfx_info = dict(zip(keys,output))
        pfx_info['mini-pfx']='{sample_id}'.format(**pfx_info)
        pfx_info['pfx']='{sample_id}'.format(**pfx_info)
        if re.search('msi',pfx_info['library-version'], re.IGNORECASE):
            pfx_info['assay']=ASSAYS[pfx_info['library-version'].upper()]
        else:
            pfx_info['assay']=ASSAYS[pfx_info['library-version']]
    elif len(output)==1:
        #Only the old LMG/OPX should hit this. 
        pfx=output[0]
        pfx_info={'mini-pfx':pfx,
                  'pfx':pfx,
                  'sample_id':pfx}
        if re.search('LMG', pfx):
            pfx_info['assay']='coloseq'
        elif re.search('OPX', pfx):
            pfx_info['assay']='oncoplex'
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
        d = datetime.datetime.strptime(date, '%y%m%d')
        return d.strftime('%Y-%m-%d')
    else:
        return date

def munge_path(pth):
    """
    Get date, run, project, machine, assay, prep-type from path
    """
    if re.search('npm1', pth.lower()):
        raise ValueError('npm1 assay not included in munge utils')
    elif re.search('msiplus', pth.lower()):
        raise ValueError('msiplus assay not included in munge utils')

    output=multi_split(pth, '/_')

    #Assuming we want YYMMDD_RUN_PROJECT
    if output[-1]=='output':
        output=output[-4:-1]
    #If old version of data that isn't in a 'output' subfolder
    else:
        output=output[-3:]
    run_pattern = re.compile('^[HMN].*_')
    #If this is correctly formatted run, process
    if re.match('\d{6}', output[0]) and re.match('^[HNM]',output[1]):
        keys=['date','run', 'project']
        pathinfo = dict(zip(keys,output))
        pathinfo['date']=munge_date(pathinfo['date'])
        #Lowercase project
        pathinfo['project']=pathinfo['project'].lower()
        #Set Machine
        pathinfo['machine']=MACHINES[pathinfo['run'][0]]
        #Set assay
        if re.search('colo', pathinfo['project']):
            pathinfo['assay']='coloseq'
        elif re.search('onco', pathinfo['project']):
            pathinfo['assay']='oncoplex'
        elif re.search('epi', pathinfo['project']):
            pathinfo['assay']='epiplex'
        elif re.search('imm', pathinfo['project']):
            pathinfo['assay']='immunoplex'
        elif re.search('marrow', pathinfo['project']):
            pathinfo['assay']='marrowseq'
        #Set prep type
        if re.search('kapa', pathinfo['project']):
            pathinfo['prep_type']='kapa'
        else:
            pathinfo['prep_type']='sure_select'
        return pathinfo
    else:
        raise ValueError('Run folder name not properly formatted')

    #def get_info(fname, pfx=None, run=None, project=None):

def munge_samples(pth):
    """
    Get pfx, run, and project from either the manifest or the sample name
    """
    #Strip trailing slash if present
    pth = os.path.normpath(pth)
    manifest=os.path.join(pth,'configs/manifest.csv')
    project=os.path.basename(pth)
    #os.path.isfile returns True if file exists
    info={}
    if os.path.isfile(manifest):        
        reader=csv.DictReader(open(manifest))
        for row in reader:
            pfx="-".join([row['run_number'], row['barcode_id']])
            is_control = True if row['sample_type'].lower()=='c' else False
            info=dict(pfx=pfx,run=row['run_number'],project=project,is_control=is_control)
            yield info
    else:
        try:
            pfx_info = munge_pfx(pth.fname)
            with open(os.path.join(pth.dir, pth.fname)) as fname:
                log.debug('loading variants from file %s' % fname)
                db._add_variants(get_info(fname, pfx=pfx_info['pfx'], run=pfx_info['machine-run'], project=run_info['project']), allow_missing=True)
        except ValueError:
            pfx = fix_pfx(pth.fname.split('_')[0])
            run = munge_path(args.path)['run']
            project = munge_path(args.path)['project'].lower()
            with open(os.path.join(args.path, pth.fname)) as fname:
                # Insert a row of data, use get_info to get the data. Repeat for all
                db._add_variants(get_info(fname, pfx=pfx, run=run, project=project), allow_missing=True)

        print "no manifest munging for you"

