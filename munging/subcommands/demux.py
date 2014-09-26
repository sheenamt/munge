"""
Crawl analysis files to create one analysis file with all info

Usage:

 munge demux flowcell-ID cores sample-sheet {narwhal,larf}

"""

# import argparse
import os
# import subprocess

# parser = argparse.ArgumentParser()
def build_parser(parser):
    parser.add_argument('flowcell-ID',
                    help="Example: 130411_D00180_0059_AH0E3KADXX")
    parser.add_argument('cores',
                    help="Number of cores for demux")
    parser.add_argument('sample-sheet',
                    help="Absolute path to sample sheet")
    parser.add_argument('server',choices=['narwhal','larf'],
                    help="Server name")

SEQ_MACHINES={'D00180':{'machine':'HA',
                        'server':'narwhal',
                        'drive':'/home/illumina/hiseq'}}

def parse_flowcell(flowcell_dir):
    """Parse the date, machine, run, side and flowcell ID
    from the flowcell directory exp:140902_D00180_0185_AHAC3AADXX"""
    fc_dict=['flowcell_dir','run_date','machine_id', 'machine_run', 'machine_side', 'flowcell_id']
    fc_values=flowcell_dir.split('_')
    fc_values.insert(0,flowcell_dir)
    fc=fc_values.pop()
    #Parse the side [A/B] from the start of the flowcell ID [A]HAC3AADXX
    fc_values.append(fc[0])
    fc_values.append(fc[1:])
    run_info=zip(fc_dict, fc_values)
    return run_info

def run_bcl2fastq(run_info, cores):
    """Run bcl2fastq for de-multiplexing and fastq generation.
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

def cat_fastqs(run_info):
    """Concatenate fastqs and write to output directory"""
    # #Now we concatenate all the fastqs together. Change "Project_default" if your project is named in the sample sheet.
    project_dirs = os.listdir(run_info['fastq_output_dir'])
    for project in project_dirs:
        if project.startswith("Project"):
            output_dir=os.path.join(run_info['drive'],run_info['fastq_output_dir']+"_"+project.split('_')[-1])
            project_dir=os.path.join(run_info['fastq_output_dir'],project)
            p1 = subprocess.Popen(["ls" ,project_dir], stdout=subprocess.PIPE) #Set up the ls command and direct the output to a pipe
            p2 = subprocess.Popen(["xargs", "-P10", "-I", "file", "cat_fastqs.sh", "file", project_dir, output_dir], stdin=p1.stdout) #send p1's output to p2
            p1.stdout.close() #make sure we close the output so p2 doesn't hang waiting for more input
            p2.communicate() #run

#args = parser.parse_args()
def action(args):
    info=vars(args)
    run_info={}
    # Parse the flowcell dir to create the run_info
    run_info.update(parse_flowcell(info['flowcell-ID']))
    # Add the sample sheet to the run dict
    run_info.update({'SampleSheet': info['sample-sheet']})
    # Add server based on seq machine id
    run_info.update(SEQ_MACHINES[run_info['machine_id']])
    run_bcl2fastq(run_info, info['cores'])
    print "run info:", run_info
    cat_fastqs(run_info)


