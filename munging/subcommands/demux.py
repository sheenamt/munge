"""
Crawl analysis files to create one analysis file with all info

Usage:

 munge demux run-folder cores -s sample-sheet 

"""

import glob
import re
import shutil
import os
import sys
import subprocess

# parser = argparse.ArgumentParser()
def build_parser(parser):
    parser.add_argument('run_folder',
                    help="Example: /home/illumina/hiseq/130411_D00180_0059_AH0E3KADXX")
    parser.add_argument('cores',
                    help="Number of cores for demux")
    parser.add_argument('-s','--samplesheet',
                    help="Absolute path to sample sheet")
    parser.add_argument('-f','--fileshare',
                    default='/mnt/disk3/genetics/fastqs',
                    help="Absolute path to fastq output on file share")

SEQ_MACHINES={'D00180':{'machine':'HA',
                        'server':'narwhal'},
              'NS500359':{'machine':'NA',
                        'server':'narwhal'},
              'M00829':{'machine':'MA',
                        'server':'narwhal'}}

def parse_flowcell_dir(run_folder):
    """Parse the date, machine, run, side and flowcell ID
    from the flowcell directory exp:
    HiSeq:   140902_D00180_0185_AHAC3AADXX
    NextSeq: 140911_NS500359_0002_AH12C0BGXX
    MiSeq:   141023_M00829_0003_000000000-AAC7L"""
    fc_dict=['flowcell_dir','run_date','Illumina_id','machine_run','machine_side','flowcell_id','drive']
    #Parse the flowcell from the run_folder
    flowcell_dir=os.path.basename(run_folder)
    #Parse the drive from the run_folder
    drive=os.path.dirname(run_folder)
    fc_values=flowcell_dir.split('_')
    fc_values.insert(0,flowcell_dir)
    fc=fc_values.pop()
    #Parse the side [A/B] from the start of the flowcell ID [A]HAC3AADXX
    if fc[0] is 'A' or fc[0] is 'B':
        fc_values.append(fc[0])
        fc_values.append(fc[1:])
    #MiSeq has a different flowcell ID format, because Illumina. 
    else:
        fc_values.append('-')
        fc_values.append(fc.split('-')[1])
    fc_values.append(drive)
    run_info=zip(fc_dict, fc_values)
    return run_info

def run_bcl2fastqv2(run_info, cores):
    """Run bcl2fastqv2 for de-multiplexing and fastq 
    generation of reads from the HiSeq, MiSeq, and NextSeq. 
    run_folder -- directory of Illumina outputs
    ss_csv -- Samplesheet CSV file describing samples.
    """
    runfolder_dir = os.path.join(run_info['drive'],run_info['flowcell_dir'])
    fastq_output_dir = os.path.join(run_info['drive'],run_info['run_date']+"_"+run_info['machine']+run_info['machine_run'])
    run_info.update({'fastq_output_dir':fastq_output_dir})
    if not os.path.exists(fastq_output_dir):
        subprocess.call(['bcl2fastq-v2.17',
                         '--runfolder-dir', runfolder_dir,
                         '--output-dir', fastq_output_dir,
                         '--sample-sheet', run_info['SampleSheet'],
                         '--processing-threads', cores,
                         '--with-failed-reads',
                         '--no-lane-splitting',
                         '--barcode-mismatches', '0',
                         '--use-bases-mask','Y*,I8,Y*'])
    return run_info

def rename_fastqs(run_info, fileshare):
    """Rename 65R21-E03-OPXv3_S8_R1_001.fastq.gz
    to 65R21_E03_OPXv3_NA0006.1.fastq.gz"""
    fastq_dirs = os.listdir(run_info['fastq_output_dir'])
    run=run_info['machine']+run_info['machine_run']
    fileshare_output= os.path.join(fileshare,run_info['run_date']+"_"+run)
    for project in fastq_dirs:
        infiles = glob.glob(os.path.join(run_info['fastq_output_dir'], project,'*fastq.gz'))
        #There are other directories in here, we only care about the ones with fastqs
        if len(infiles)>=2:
            project_dir = '_'.join([fileshare_output,project])
            if not os.path.exists(project_dir):
                os.makedirs(project_dir)
            for f in infiles:
                pfx = os.path.basename(f).split('_')[0]
                print pfx
                if re.search('_R1_', f):
                    #Move newly named fastq to DATE_RUN_PROJECT directory
                    newname = os.path.join(project_dir, '%s_%s.1.fastq.gz' % (pfx.replace('-','_'), run))
                    shutil.copy2(f, newname)
                elif re.search('_R2_', f):
                    #Move newly named fastq to DATE_RUN_PROJECT directory
                    newname = os.path.join(project_dir, '%s_%s.2.fastq.gz' % (pfx.replace('-','_'), run))
                    shutil.copy2(f, newname)


def action(args):
    info=vars(args)
    run_info={}
    # Parse the flowcell dir to create the run_info
    run_info.update(parse_flowcell_dir(info['run_folder']))
    # Add the sample sheet to the run dict
    if args.samplesheet:
        run_info.update({'SampleSheet': info['samplesheet']})
    else:
        run_info.update({'SampleSheet': os.path.join(info['run_folder'],'SampleSheet.csv')})
    # Add server based on seq machine id
    run_info.update(SEQ_MACHINES[run_info['Illumina_id']])
    run_bcl2fastqv2(run_info, info['cores'])
    print "run info:", run_info
    rename_fastqs(run_info, args.fileshare)

    print "ran /home/local/AMC/sheenams/bcl2fastq"
