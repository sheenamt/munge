"""
Run bcl2fastq to demultiplex a sequencing run

Usage:

 munge demux run-folder cores sequencer

"""

# import argparse
import os
import sys
import subprocess
import glob

# parser = argparse.ArgumentParser()
def build_parser(parser):
    parser.add_argument('run_folder',
                    help="Example: /home/illumina/hiseq/130411_D00180_0059_AH0E3KADXX")
    parser.add_argument('cores',
                    help="Number of cores for demux")
    parser.add_argument('-s','--samplesheet',
                    help="Absolute path to sample sheet")
    parser.add_argument('sequencer',choices=['HiSeq','MiSeq','NextSeq'],
                    help="Sequencer name")

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
    generation of reads from the HiSeq, MiSeq, and NextSeq. 
    run_folder -- directory of Illumina outputs
    ss_csv -- Samplesheet CSV file describing samples.
    """
    runfolder_dir = os.path.join(run_info['drive'],run_info['flowcell_dir'])
    fastq_output_dir = os.path.join(run_info['drive'],run_info['run_date']+"_"+run_info['machine']+run_info['machine_run'])
    run_info.update({'fastq_output_dir':fastq_output_dir})
    if not os.path.exists(fastq_output_dir):
        subprocess.call(['/home/genetics/working/sheenams/demux-v2/tmp/bcl2fastq2-v2.16.0.10-build/bcl2fastq2-v2.16.0.10/bin/bcl2fastq',
                         '--runfolder-dir', runfolder_dir,
                         '--output-dir', fastq_output_dir,
                         '--sample-sheet', run_info['SampleSheet'],
                         '--processing-threads', cores,
                         '--with-failed-reads',
                         '--barcode-mismatches', '0',
                         '--no-lane-splitting',
                         '--use-bases-mask','Y*,I8,Y*'])

    return run_info

# def run_bcl2fastqv2(run_info, cores):
#     """Run bcl2fastqv2 for de-multiplexing and fastq 
#     generation of reads from the NextSeq. 
#     run_folder -- directory of Illumina outputs
#     ss_csv -- Samplesheet CSV file describing samples.
#     """
#     runfolder_dir = os.path.join(run_info['drive'],run_info['flowcell_dir'])
#     fastq_output_dir = os.path.join(run_info['drive'],run_info['run_date']+"_"+run_info['machine']+run_info['machine_run'])
#     print fastq_output_dir
#     run_info.update({'fastq_output_dir':fastq_output_dir})
#     if not os.path.exists(fastq_output_dir):
#         subprocess.call(['bcl2fastq',
#                          '--runfolder-dir', runfolder_dir,
#                          '--output-dir', fastq_output_dir,
#                          '--demultiplexing-threads', cores,
#                          '--processing-threads', cores,
#                          '--with-failed-reads',
#                          '--use-bases-mask','Y*,I8,Y*'])
#     return run_info

def cat_fastqs(run_info):
    """Concatenate fastqs and write to output directory"""
    # #Now we concatenate all the fastqs together. Change "Project_default" if your project is named in the sample sheet.
    project_dirs = os.listdir(run_info['fastq_output_dir'])
    for project in project_dirs:
        if project.startswith("Project"):
            print "project:", project
            output_dir=os.path.join(run_info['drive'],run_info['fastq_output_dir']+"_"+project.split('_')[-1])
            print "output_dir:", output_dir
            project_dir=os.path.join(run_info['fastq_output_dir'],project)
            print "project_dir:", project_dir
            seq_run = ''.join([run_info['machine'],run_info['machine_run']])
            files =glob.glob(os.path.join(run_info['fastq_output_dir'],project, "*_R1_*"))
            for f in files:
                print f
                p2 = subprocess.Popen(["cat_fastqs_v2.sh", f, project_dir, output_dir, seq_run] ) #stdin=p1.stdout) #send p1's output to p2
#            p1 = subprocess.Popen(['ls *_R1_*', 'project_dir']) #, stdout=subprocess.PIPE) 
#            p1.communicate()
            #Set up the ls command and direct the output to a pipe

#            p1.stdout.close() #make sure we close the output so p2 doesn't hang waiting for more input
            p2.communicate() #run

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
    # if info['sequencer'] == 'NextSeq':
    #     print "running bcl2fastq on nextseq data"
    #     run_bcl2fastqv1(run_info, info['cores'])
    # else:


    #concatenate the fastqs across lanes
    cat_fastqs(run_info)
    print "run info:", run_info
