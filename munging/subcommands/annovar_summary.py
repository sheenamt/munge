"""
Summarize output from Annovar and EVS

Usage:

 munge summary /path/to/captured/genes/ $SAVEPATH/$PFX.* -o $SAVEPATH/${PFX}_Analysis.txt;

"""

import logging
import sys
import re
import argparse
import csv
from collections import defaultdict
from os import path

from munging.annotation import get_location, multi_split, split_string_in_two

variant_headers = ['chr', 'start', 'stop', 'Ref_Base', 'Var_Base']

# (pattern, header_ids, var_key_ids)
file_types = {
#files
    'variant_function': ({0: 'var_type_1',
                          1: 'Gene',
                          7: 'Zygosity',
                          12: 'rsid_1',
                          18: 'Read_Headers',
                          19: 'Reads',
                          8: 'GATK_Score'},
                         [2, 3, 4, 5, 6]),
    'exonic_variant_function': ({1: 'var_type_2', 2: 'Transcripts'}, [3, 4, 5, 6, 7]),
    'hg19_ALL.sites.2015_08_dropped': ({1: '1000g_ALL'}, [2, 3, 4, 5, 6]),
    'hg19_AMR.sites.2015_08_dropped': ({1: '1000g_AMR'}, [2, 3, 4, 5, 6]),
    'hg19_AFR.sites.2015_08_dropped': ({1: '1000g_AFR'}, [2, 3, 4, 5, 6]),
    'hg19_SAS.sites.2015_08_dropped': ({1: '1000g_SAS'}, [2, 3, 4, 5, 6]),
    'hg19_EAS.sites.2015_08_dropped': ({1: '1000g_EAS'}, [2, 3, 4, 5, 6]),
    'hg19_EUR.sites.2015_08_dropped': ({1: '1000g_EUR'}, [2, 3, 4, 5, 6]),
    'hg19_exac03_dropped': ({1: 'EXAC'}, [2, 3, 4, 5, 6]),
    'hg19_cosmic70_dropped': ({1: 'Cosmic'}, [2, 3, 4, 5, 6]),
    'hg19_genomicSuperDups': ({0: 'Segdup'}, [2, 3, 4, 5, 6]),
    'hg19_dbnsfp30a_dropped': ({1: 'ljb_Scores'}, [2, 3, 4, 5, 6]),
    'hg19_esp6500siv2_all_dropped': ({1: 'EVS_esp6500_ALL'}, [2, 3, 4, 5, 6]),
    'hg19_esp6500siv2_ea_dropped': ({1: 'EVS_esp6500_EU'}, [2, 3, 4, 5, 6]),
    'hg19_esp6500siv2_aa_dropped': ({1: 'EVS_esp6500_AA'}, [2, 3, 4, 5, 6]),
    'hg19_UW_freq_dropped': ({1: 'UW_Freq_list'}, [2, 3, 4, 5, 6]),
    'hg19_variants_dropped': ({1: 'Clinically_Flagged'}, [2, 3, 4, 5, 6]),
    'hg19_nci60_dropped': ({1: 'NCI60'}, [2, 3, 4, 5, 6]),
    'hg19_clinvar_20150629_dropped': ({1: 'ClinVar'}, [2, 3, 4, 5, 6]),
    'hg19_CADD_dropped': ({1: 'CADD'}, [2, 3, 4, 5, 6]),
    'hg19_snp138_dropped': ({1: 'rsid_2'}, [2, 3, 4, 5, 6]),
    'hg19_dbscsnv11_dropped': ({1: 'splicing'}, [2, 3, 4, 5, 6]), #probability score for each variant that reflects the confidence that the variant alters splicing.
    'hg19_clinical_variants_dropped': ({1: 'Clinically_Flagged'}, [2, 3, 4, 5, 6]),
}
log = logging.getLogger(__name__)


def get_reads(headers, data):
    """Parse the reads from
    GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR
    OR
    GT:AD:DP:GQ:PL
    RD:Depth of reference-supporting bases (reads1)
    AD:Depth of variant-supporting bases (reads2) OR (reads1,reads2)
    ABQ:Average quality of variant-supporting bases (qual2)
    """
    info = dict(zip(headers.split(':'), data.split(':')))
    #Do not return GATK reads. They are downsampled to 250 
    if 'RD' in info.keys():
        return info['RD'],info['AD'],info['ABQ']
    else:
        return '-1','-1',''

        


def munge_gene_and_Transcripts(data, RefSeqs):
    """
    Return modified values of (Gene, Transcripts). Note that
    this depends on 'Variant_Type' provided by munge_variant.
    """
    Transcripts = data.get('Transcripts')
    Gene = data.get('Gene', '')
    if not Gene or data['Variant_Type'] in ('upstream','downstream','intergenic','ncRNA_exonic'):
        Gene = ''
    elif '(' in Gene:
        Gene, Gene_tail = data['Gene'].split('(', 1)
        if 'NM' in Gene_tail:
            # overrides value of Transcripts
            Transcripts = data['Gene'].replace('(', ':').strip(')')

    return Gene, Transcripts


def munge_transcript(data, RefSeqs):
    """
    Return HGVS correct transcript annotations
    Filtered with a preferred transcript list
    NM_006772.1:c.1713G>A
    """
    codon, prot, protein, coding, txpt = ' ', ' ', ' ', ' ', None
    transcripts = data.get('Transcripts')
    if transcripts is not None:
        # Split incoming trans, strip the trailing )
        data = multi_split(transcripts.replace('(', ':'), ',(')
        for d in data:
            # Split the actual transcript info which is colon separated
            x = d.split(':')
            try:
                gene, txpt, exon, codon, prot = x

            except ValueError:
                try:
                    gene, txpt, exon, codon = x
                except ValueError:
                    continue
            pref_trans = RefSeqs.get(txpt)
            #Want to return None for all values if not pref_trans
            if not pref_trans:
                continue
            code = pref_trans + ':' + codon
            coding = coding + ' ' + code
            protein = protein + ' ' + prot
    return coding.strip(), protein.strip()


def map_headers(fname, header_ids, variant_idx):
    """
    Gets header(s) and info from each file
    """
    # dict mapping variant keys to columns
    with open(fname, 'rU') as infile:
        reader = csv.reader(infile, delimiter='\t')
        for row in reader:

            #adds the value from each position (2,3,4,5,6) to the variant_id
            variant_id = tuple(row[i] for i in variant_idx)

            #create dict of variants {2:'chr', 3:'start', 4:'stop, 5:'ref', 6:'var'}
            variant_dict = dict(zip(variant_headers, variant_id))

            #create dict of data from this file; key is included only
            #if value is not the empty string
            data = dict((key, row[i]) for i, key in header_ids.items() if row[i].strip())
            data['Position'] = get_location(**variant_dict)

            yield (variant_id, dict(data, **variant_dict))


def get_allele_freq(data):
    """
    Return allele frequency of var_reads/ref_reads
    """
    if int(data['Var_Reads']) == -1 and int(data['Ref_Reads']) == -1:
        freq = 'NA'
    else:
        total = int(data['Var_Reads']) + int(data['Ref_Reads'])
        freq = int(data['Var_Reads']) / float(total)
        freq = "{0:.2f}".format(freq)
    return freq

def munge_ljb_scores(data):
    """
    Parse sift, polyphen and gerp from ljb_all file
    """
    try:
        info = data.get('ljb_Scores').split(',')
        polyphen = info[2]
        sift = info[0]
        gerp = info[28]
        mutation_taster = info[8]
    except AttributeError:
        return -1, -1, -1, -1 
    return polyphen, sift, mutation_taster, gerp

def build_parser(parser):
    parser.add_argument('RefSeqs', type=argparse.FileType('rU'),
        help='Capture genes file with RefSeq in second column')
    parser.add_argument(
        'infiles', action='append', nargs='+',
        help='Input files')
    parser.add_argument(
        '-o', '--outfile',
        help='Output file', default=sys.stdout,
        type=argparse.FileType('w'))
    parser.add_argument(
        '--strict', action='store_true', default=False,
        help='Exit with error if an input file has no match.')


def action(args):

    (infiles, ) = args.infiles

    RefSeqs = {}
    if args.RefSeqs:
        refs = csv.DictReader(args.RefSeqs, delimiter='\t')
        for row in refs:
            if row['RefSeq'] :
                for transcript in row['RefSeq'].split('/'):
                    RefSeqs[transcript.split('.')[0]] = transcript

    headers = ['Position'] + variant_headers[3:5] + [
        'Clinically_Flagged',
        'Variant_Type',
        'UW_Freq',
        'Filter',
        '1000g_ALL',
        'EVS_esp6500_ALL',        
        'EXAC',
        'Gene',
        'p.',
        'c.',
        'Faves_Y/N',
        'Ref_Reads',
        'Var_Reads',
        'Allele_Frac',
        'Variant_Phred',
        'Cosmic',
        'CADD',
        'ClinVar',
        'Polyphen',
        'Sift',
        'Mutation_Taster',
        'Gerp',
        '1000g_AMR',
        '1000g_EUR',
        '1000g_SAS',
        '1000g_EAS',
        '1000g_AFR',
        'EVS_esp6500_AA',
        'EVS_esp6500_EU',
        'Transcripts',
        'Zygosity',
        'Segdup',
        'NCI60',
        'dbSNP_ID',
        'UW_Count',
        'GATK_Score',
        'ADA_Alter_Splice',
        'RF_Alter_Splice',
        ]

    writer = csv.DictWriter(args.outfile,
                            fieldnames=headers,
                            quoting=csv.QUOTE_MINIMAL,
                            extrasaction='ignore',
                            delimiter='\t')

    writer.writeheader()

    # accumulate data from all input files for each variant
    output = defaultdict(dict)
    for fname in infiles:
        try:
            _, file_type = path.basename(fname).split('.', 1)
        except ValueError:
            continue
        try:
            #if the file type matches one we want,
            #header ids are output columns
            #var_key_ids are chrm:str:stp:ref:var
            header_ids, var_key_ids = file_types[file_type]
        except KeyError:
            if re.search('dropped', fname):
                log.warning('no match: %s' % fname)
            if args.strict:
                sys.exit(1)
            continue
        multi_trans_keys=['Transcripts','var_type_1','var_type_2']
        for var_key, data in map_headers(fname, header_ids, var_key_ids):
            #If position already in output{}, update certain fields
            if var_key in output:
                for key in multi_trans_keys:
                    if key in output[var_key].keys() and data.get(key):
                        if data.get(key) not in output[var_key][key]:
                            output[var_key][key]=','.join([output[var_key][key],data.get(key)])
                    elif key in output[var_key].keys() and not data.get(key):
                        output[var_key][key]=output[var_key][key]
                    elif key not in output[var_key].keys() and data.get(key):
                        output[var_key][key]=data.get(key)
                #grab the keys of data, removing keys list, and update those keys
                data_keys=data.keys()
                if output[var_key].has_key('Reads'):
                    if 'RD' in output[var_key]['Read_Headers']:
                        output[var_key]['Read_Headers']=output[var_key]['Read_Headers']
                        output[var_key]['Reads']=output[var_key]['Reads']
                        data_keys.remove('Reads')
                        data_keys.remove('Read_Headers')
                    else:
                        output[var_key]['Read_Headers']=data.get('Read_Headers')
                        output[var_key]['Reads']=data.get('Reads')
                  
                for k in data_keys:
                    if k not in multi_trans_keys:
                        output[var_key][k]=data.get(k)
            else:
                output[var_key].update(data)
    sort_key = lambda row: [(row[k]) for k in ['chr', 'start', 'stop', 'Ref_Base', 'Var_Base']]
    # # write each row (with all data aggregated), modifying fields as necessary
    for data in sorted(output.values(), key=sort_key):
        # # # modify any specific fields here
        data['Variant_Type'] = data.get('var_type_2') if data.get('var_type_2') else data.get('var_type_1')
        data['Gene'], data['Transcripts'] = munge_gene_and_Transcripts(data, RefSeqs)
        data['c.'], data['p.'] = munge_transcript(data, RefSeqs)
        data['Polyphen'], data['Sift'],data['Mutation_Taster'],data['Gerp'] = munge_ljb_scores(data)
        data['dbSNP_ID'] = data.get('rsid_1') or data.get('rsid_2')
        data['1000g_ALL'] = data.get('1000g_ALL') or -1
        data['1000g_AMR'] = data.get('1000g_AMR') or -1
        data['1000g_SAS'] = data.get('1000g_SAS') or -1
        data['1000g_EAS'] = data.get('1000g_EAS') or -1
        data['1000g_AFR'] = data.get('1000g_AFR') or -1
        data['1000g_EUR'] = data.get('1000g_EUR') or -1
        data['EXAC'] = data.get('EXAC').split(',')[0] if data.get('EXAC') else -1      
        data['EVS_esp6500_ALL'] = data.get('EVS_esp6500_ALL').split(',')[0] if data.get('EVS_esp6500_ALL') else -1
        data['EVS_esp6500_AA'] = data.get('EVS_esp6500_AA').split(',')[0] if data.get('EVS_esp6500_AA') else -1
        data['EVS_esp6500_EU'] = data.get('EVS_esp6500_EU').split(',')[0] if data.get('EVS_esp6500_EU') else -1
        #CADD is raw score, phred score. We only care about phred
        _, data['CADD'] = split_string_in_two(data.get('CADD'))
        data['ADA_Alter_Splice'],data['RF_Alter_Splice'] = split_string_in_two(data.get('splicing'))
        data['UW_Freq'], data['UW_Count'] = split_string_in_two(data.get('UW_Freq_list'))
        data['Ref_Reads'], data['Var_Reads'], data['Variant_Phred'] = get_reads(data.get('Read_Headers'),data.get('Reads'))
        data['Allele_Frac'] = get_allele_freq(data)
        writer.writerow(data)
