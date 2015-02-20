"""
Crawl analysis files to create one analysis file with all info

Usage:

 munge demux flowcell-ID-folder cores sample-sheet sequencer

"""

# import argparse
import os
import sys
import subprocess

# parser = argparse.ArgumentParser()
def build_parser(parser):
    parser.add_argument('flowcell-ID-folder',
                    help="Example: 130411_D00180_0059_AH0E3KADXX")
    parser.add_argument('cores',
                    help="Number of cores for demux")
    parser.add_argument('sample-sheet',
                    help="Absolute path to sample sheet")
    parser.add_argument('sequencer',choices=['HiSeq','MiSeq','NextSeq'],
                    help="Sequencer name")

SEQ_MACHINES={'D00180':{'machine':'HA',
                        'server':'narwhal',
                        'drive':'/home/illumina/hiseq'},
              'NS500359':{'machine':'NA',
                        'server':'narwhal',
                        'drive':'/home/illumina/NextSeq'},
              'M00829':{'machine':'MA',
                        'server':'narwhal',
                        'drive':'/home/illumina/miseq'}}

def parse_flowcell_dir(flowcell_dir):
    """Parse the date, machine, run, side and flowcell ID
    from the flowcell directory exp:
    HiSeq:   140902_D00180_0185_AHAC3AADXX
    NextSeq: 140911_NS500359_0002_AH12C0BGXX
    MiSeq:   141023_M00829_0003_000000000-AAC7L"""
    fc_dict=['flowcell_dir','run_date','machine_id', 'machine_run', 'machine_side', 'flowcell_id']
    print flowcell_dir
    fc_values=flowcell_dir.split('_')
    fc_values.insert(0,flowcell_dir)
    fc=fc_values.pop()
    #Parse the side [A/B] from the start of the flowcell ID [A]HAC3AADXX
    if fc[0] is 'A' or fc[0] is 'B':
        fc_values.append(fc[0])
        fc_values.append(fc[1:])
    #MiSeq has a different flowcell ID format, because Illumina. 
    else:
        fc_values.append('NA')
        fc_values.append(fc.split('-')[1])
    run_info=zip(fc_dict, fc_values)
    return run_info

def run_bcl2fastqv1(run_info, cores):
    """Run bcl2fastq v1 for de-multiplexing and fastq 
    generation of reads from the HiSeq or MiSeq.
    run_folder -- directory of Illumina outputs
    ss_csv -- Samplesheet CSV file describing samples.
    """
    bc_dir = os.path.join(run_info['drive'],run_info['flowcell_dir'], "Data", "Intensities", "BaseCalls")
    fastq_output_dir = os.path.join(run_info['drive'],run_info['run_date']+"_"+run_info['machine']+run_info['machine_run'])
    run_info.update({'fastq_output_dir':fastq_output_dir})
    if not os.path.exists(fastq_output_dir):
        subprocess.call(['configureBclToFastq.pl',
                         '--input-dir', bc_dir,
                         '--output-dir', fastq_output_dir,
                         '--sample-sheet', run_info['SampleSheet'],
                         '--flowcell-id', run_info['flowcell_id'],
                         '--no-eamss',
                         '--with-failed-reads',
                         '--use-bases-mask','Y*,I8,Y*'])

        os.chdir(fastq_output_dir)
        subprocess.call(["make", "-j", cores])

    return run_info

def run_bcl2fastqv2(run_info, cores):
    """Run bcl2fastqv2 for de-multiplexing and fastq 
    generation of reads from the NextSeq. 
    run_folder -- directory of Illumina outputs
    ss_csv -- Samplesheet CSV file describing samples.
    """
    runfolder_dir = os.path.join(run_info['drive'],run_info['flowcell_dir'])
    fastq_output_dir = os.path.join(run_info['drive'],run_info['run_date']+"_"+run_info['machine']+run_info['machine_run'])
    print fastq_output_dir
    run_info.update({'fastq_output_dir':fastq_output_dir})
    if not os.path.exists(fastq_output_dir):
        subprocess.call(['bcl2fastq',
                         '--runfolder-dir', runfolder_dir,
                         '--output-dir', fastq_output_dir,
                         '--demultiplexing-threads', cores,
                         '--processing-threads', cores,
                         '--with-failed-reads',
                         '--use-bases-mask','Y*,I8,Y*'])
    return run_info

def combine_fastq_files(in_files, project_dir, output_dir):

    print in_files, project_dir, output_dir
    print os.path.basename(in_files[0]).strip('Sample_')
    for sample in in_files:
        pfx=sample.strip('Sample_')
        for R in (1,2):
            to_zip="/".join([project_dir,sample])
            zip_to="/".join([output_dir,pfx])
            p1 = subprocess.Popen(["zcat" ,to_zip,"*$R*.fastq.gz"], stdout=subprocess.PIPE) #Set up the ls command and direct the output to a pipe
#            p2 = subprocess.Popen(["|","gzip", "-n>", zip_to,"*R.fastq.gz"], stdin=p1.stdout) #send p1's output to p2

def cat_fastqs(run_info):
    """Concatenate fastqs and write to output directory"""
    # #Now we concatenate all the fastqs together. Change "Project_default" if your project is named in the sample sheet.
    project_dirs = os.listdir(run_info['fastq_output_dir'])
    for project in project_dirs:
        print "project:", project
        if project.startswith("Project"):
            output_dir=os.path.join(run_info['drive'],run_info['fastq_output_dir']+"_"+project.split('_')[-1])
            project_dir=os.path.join(run_info['fastq_output_dir'],project)
            seq_run = ''.join([run_info['machine'],run_info['machine_run']])
            p1 = subprocess.Popen(["ls" ,project_dir], stdout=subprocess.PIPE) #Set up the ls command and direct the output to a pipe
            p2 = subprocess.Popen(["xargs", "-P10", "-I", "file", "cat_fastqs.sh", "file", project_dir, output_dir, seq_run], stdin=p1.stdout) #send p1's output to p2
            p1.stdout.close() #make sure we close the output so p2 doesn't hang waiting for more input
            p2.communicate() #run

def action(args):
    info=vars(args)
    run_info={}
    # Parse the flowcell dir to create the run_info
    run_info.update(parse_flowcell_dir(info['flowcell-ID-folder']))
    # Add the sample sheet to the run dict
    run_info.update({'SampleSheet': info['sample-sheet']})
    # Add server based on seq machine id
    run_info.update(SEQ_MACHINES[run_info['machine_id']])
    if info['sequencer'] == 'NextSeq':
        run_bcl2fastqv2(run_info, info['cores'])
    else:
        run_bcl2fastqv1(run_info, info['cores'])
    print "run info:", run_info
    cat_fastqs(run_info)
