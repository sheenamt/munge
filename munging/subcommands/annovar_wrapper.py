"""
Run annovar to generate standard set of annotations to bring into the DB using annotation_importer
"""
import logging
from collections import namedtuple
import sys
import os
import subprocess

def build_parser(parser):
    parser.add_argument('out_dir',
                        help='Directory where output files will be created')
    parser.add_argument('input_file', default=None,
                        help='Explicitly specify input file of variants in Annovar format')
    parser.add_argument('--machine',
                        choices=['H', 'N', 'M'],
                        help='Machine type for internal frequencies')
    parser.add_argument('--assay',
                        help='Assay code for interal frequencies')
    parser.add_argument('--ref_gene_only', action='store_true',
                        help='Run only the ref_gene annotation')
    parser.add_argument('--clin_flagged_only', action='store_true',
                        help='Run only the clinical variants annotation')
    parser.add_argument('--clinically_flagged', default='/mnt/disk2/com/Genomes/Annovar_files/hg19_clinical_variants',
                        help='Clinically flagged variants file')
    parser.add_argument('--library_dir', default='/mnt/disk2/com/Genomes/Annovar_files',
                        help='Directory holding Annovar library files')
    parser.add_argument('--annovar_bin', default='',
                        help='Location of the Annovar perl executables')

log = logging.getLogger(__name__)


ANNOTATIONS = [('snp138',),  # dbsnp
               ('exac03',),  # ExAC 65000 exome allele frequency data 
               ('dbscsnv11',), # dbscSNV version 1.1 for splice site prediction by AdaBoost and Random Forest
               ('1000g2015aug_all',),  # 1000 genomes annotation:
               ('1000g2015aug_amr',),  # 1000 genomes (admmixed american) annotation:
               ('1000g2015aug_eur',),  # 1000 genomes (european) annotation:
               ('1000g2015aug_eas',),  # 1000 genomes (east asian) annotation:
               ('1000g2015aug_sas',),  # 1000 genomes (south asian) annotation:
               ('1000g2015aug_afr',),  # 1000 genomes (african) annotation:
               ('dbnsfp33a',),  # whole-exome SIFT, PolyPhen2 HDIV, PolyPhen2 HVAR, LRT, MutationTaster, MutationAssessor, FATHMM, MetaSVM, MetaLR, VEST, CADD, GERP++, PhyloP and SiPhy scores from dbNSFP version 3.3
               ('esp6500siv2_all',),  # alternative allele frequency in the NHLBI-ESP project with 6500 exomes, including the indel calls and the chrY calls. evs-all
               ('esp6500siv2_aa',),  # evs-african american
               ('esp6500siv2_ea',),  # evs-european
               ('cosmic84',),  # cosmic84, created internally
               ('clinvar_20180603',),  # CLINVAR database with Variant Clinical Significance (unknown, untested, non-pathogenic, probable-non-pathogenic, probable-pathogenic, pathogenic, drug-response, histocompatibility, other) and Variant disease name, Clinvar version 20170905 with separate columns (CLINSIG CLNDBN CLNACC CLNDSDB CLNDSDBID
               ('nci60',),  # NCI-60 human tumor cell line panel exome sequencing allele frequency data
               ('segdup', '--regionanno',),  # segdup annotation:              
               ('refGene', '--geneanno', ['--splicing_threshold','10', '--hgvs']),  # Gene level annotation:
#               ('intervar_20180118',), #clinical interpretation of missense variants (indels not supported)
#               ('cg69',), #allele frequency in 69 human subjects sequenced by Complete Genomics
#               ('gnomad_exome',), #gnomAD exome collection
#               ('gnomad_genome',), #gnomAD genome collection
#               ('kaviar_20150923',), #170 million Known VARiants from 13K genomes and 64K exomes in 34 projects
#               ('hrcr1',), #40 million variants from 32K samples in haplotype reference consortium
#               ('abraom',), #2.3 million Brazilian genomic variants
#               ('gme',), #Great Middle East allele frequency including NWA (northwest Africa), NEA (northeast Africa), AP (Arabian peninsula), Israel, SD (Syrian desert), TP (Turkish peninsula) and CA (Central Asia)
#               ('mcap',), #M-CAP scores for non-synonymous variants
#               ('revel',), #REVEL scores for non-synonymous variants
#               ('avsnp150',), #dbSNP150 with allelic splitting and left-normalization
#               ('icgc21',), # International Cancer Genome Consortium version 21
#               ('popfreq_all_20150413',), #A database containing all allele frequency from 1000G, ESP6500, ExAC and CG46
#               ('mitimpact24',), # pathogenicity predictions of human mitochondrial missense variants
#               ('gerp++elem',), #conserved genomic regions by GERP++
#               ('cadd13',), #CADD version 1.3
#               ('fathmm',), #whole-genome FATHMM_coding and FATHMM_noncoding scores
#               ('gwava',), # whole genome GWAVA_region_score and GWAVA_tss_score
#               ('eigen',), # whole-genome Eigen scores
           ]
        

def action(args):
    BUILDVER = 'hg19'
    ANNOVAR_VARIANTS = os.path.join(args.annovar_bin, 'annotate_variation.pl')

    # Named tuple for parsing annotation defs
    AnnotInfo =  namedtuple('AnnotInfo', ['dbtype', 'anno_type', 'args', 'library_dir'])
    AnnotInfo.__new__.__defaults__ = ('-filter','', args.library_dir)

    variants_file = args.input_file
    file_pfx = os.path.basename(args.input_file).replace('.ann','')

    #Has default file, so will always be set
    generic_db_fname = os.path.basename(args.clinically_flagged)
    generic_db_dirname = os.path.dirname(args.clinically_flagged)

    #Sometimes all we want is the variant_function output 
    if args.ref_gene_only:
        annots = [AnnotInfo(dbtype='refGene', anno_type='--geneanno', args=['--splicing_threshold', '10', '--hgvs'], library_dir=args.library_dir)]
    elif args.clin_flagged_only:
        annots = [AnnotInfo(dbtype='generic', anno_type='--genericdbfile', args=[generic_db_fname, '-filter'], library_dir=generic_db_dirname)]
    else:
        annots = [AnnotInfo(*a) for a in ANNOTATIONS]
        #Machine and assay are optional and only used when we want full annotation 
        if args.machine:
            machine=args.machine
        if args.assay:
            assay=args.assay.split('v')[0]
            internal_cadd_file = '_'.join(['hg19','CADD',assay])
            if os.path.isfile(os.path.join(args.library_dir, internal_cadd_file)):
                annots.append(AnnotInfo(dbtype='generic', anno_type='--genericdbfile', args=[internal_cadd_file, '-filter']))

        if args.machine and args.assay:
            internal_freq_file = '_'.join(['hg19',machine,assay])
            if os.path.isfile(os.path.join(args.library_dir, internal_freq_file)):
                annots.append(AnnotInfo(dbtype='generic', anno_type='--genericdbfile', args=[internal_freq_file, '-filter']))
            
        #There is always a clin flagged
        annots.append(AnnotInfo(dbtype='generic', anno_type='--genericdbfile', args=[generic_db_fname, '-filter'], library_dir=generic_db_dirname))

    # run annotate_variants based on the spec in ANNOTATIONS
    for a in annots:
        annovar_cmd = [ANNOVAR_VARIANTS, 
                       '-dbtype', a.dbtype,
                       '--buildver', BUILDVER, 
                       a.anno_type ] \
                       + list(a.args) + \
                       ['--otherinfo', 
                        '--separate', 
                        '-outfile', os.path.join(os.path.dirname(args.out_dir),file_pfx),
                        variants_file, 
                        a.library_dir]
        cmd = filter(None, annovar_cmd)

        subprocess.check_call(cmd)
        if a.anno_type=='--genericdbfile':
            generic_file = os.path.join(os.path.dirname(args.out_dir),file_pfx+'.hg19_generic_dropped')
            #Case for moving frequency file
            if args.machine and (args.machine in a.args[0]):
                gen_file_basename = a.args[0].replace(machine+'_'+assay,'UW_freq')
            #Case for moving CADD file and cli
            elif args.assay and (args.assay in a.args[0]):
                gen_file_basename = a.args[0].replace(assay,'').strip('_')
            else:
                gen_file_basename = generic_db_fname
            specific_file = os.path.join(os.path.dirname(args.out_dir),file_pfx+'.'+gen_file_basename+'_dropped')
            mvcmd=['mv' , generic_file, specific_file]
            subprocess.check_call(mvcmd)




