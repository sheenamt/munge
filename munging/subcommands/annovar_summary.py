"""
Summarize output from Annovar and EVS

Usage:

 munge annovar_summary /path/to/captured/genes/ $SAVEPATH/$PFX.* -o $SAVEPATH/${PFX}_Analysis.txt;

"""

import logging
import sys
import re
import argparse
import csv
from collections import defaultdict
from os import path

from munging.annotation import get_location, multi_split, split_string_in_two



def build_parser(parser):
    parser.add_argument('RefSeqs', type=argparse.FileType('rU'),
        help='Capture genes file with RefSeq in second column')
    parser.add_argument(
        'type', choices=['SNP', 'INDEL', 'PINDEL'],
        help='Type of files to create tab, SNP or INDEL')
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

variant_headers = ['chr', 'start', 'stop', 'Ref_Base', 'Var_Base']

# (pattern, header_ids, var_key_ids)
snp_file_types = {
#files
    'merged.variant_function': ({0: 'var_type_1',
                          1: 'Gene',
                          7: 'Zygosity',
                          12: 'rsid_1',
                          18: 'Read_Headers',
                          19: 'Reads',
                          8: 'GATK_Score'},
                         [2, 3, 4, 5, 6]),
    'merged.hg19_ALL.sites.2015_08_dropped': ({1: '1000g_ALL'}, [2, 3, 4, 5, 6]),
    'merged.hg19_AMR.sites.2015_08_dropped': ({1: '1000g_AMR'}, [2, 3, 4, 5, 6]),
    'merged.hg19_AFR.sites.2015_08_dropped': ({1: '1000g_AFR'}, [2, 3, 4, 5, 6]),
    'merged.hg19_SAS.sites.2015_08_dropped': ({1: '1000g_SAS'}, [2, 3, 4, 5, 6]),
    'merged.hg19_EAS.sites.2015_08_dropped': ({1: '1000g_EAS'}, [2, 3, 4, 5, 6]),
    'merged.hg19_EUR.sites.2015_08_dropped': ({1: '1000g_EUR'}, [2, 3, 4, 5, 6]),
    'merged.hg19_exac03_dropped': ({1: 'EXAC'}, [2, 3, 4, 5, 6]),
    'merged.hg19_cosmic92_dropped': ({1: 'Cosmic'}, [2, 3, 4, 5, 6]),
    'merged.hg19_genomicSuperDups': ({0: 'Segdup'}, [2, 3, 4, 5, 6]),
    'merged.hg19_dbnsfp30a_dropped': ({1: 'ljb_Scores'}, [2, 3, 4, 5, 6]),
    'merged.hg19_esp6500siv2_all_dropped': ({1: 'EVS_esp6500_ALL'}, [2, 3, 4, 5, 6]),
    'merged.hg19_esp6500siv2_ea_dropped': ({1: 'EVS_esp6500_EU'}, [2, 3, 4, 5, 6]),
    'merged.hg19_esp6500siv2_aa_dropped': ({1: 'EVS_esp6500_AA'}, [2, 3, 4, 5, 6]),
    'merged.hg19_UW_freq_dropped': ({1: 'UW_Freq_list'}, [2, 3, 4, 5, 6]),
    'merged.hg19_nci60_dropped': ({1: 'NCI60'}, [2, 3, 4, 5, 6]),
    'merged.hg19_clinvar_20200316_dropped': ({1: 'ClinVar'}, [2, 3, 4, 5, 6]),
    'merged.hg19_CADD_dropped': ({1: 'CADD'}, [2, 3, 4, 5, 6]),
    'merged.hg19_snp138_dropped': ({1: 'rsid_2'}, [2, 3, 4, 5, 6]),
    'merged.hg19_dbscsnv11_dropped': ({1: 'splicing'}, [2, 3, 4, 5, 6]), #probability score for each variant that reflects the confidence that the variant alters splicing.
    'merged.hg19_clinical_variants_dropped': ({1: 'Clinically_Flagged'}, [2, 3, 4, 5, 6]),
    'merged.exonic_variant_function': ({1: 'var_type_2', 2: 'Transcripts'}, [3, 4, 5, 6, 7]),
    'uw_dec_p_values.csv': ({1: 'UW_DEC_p'}, [2, 3, 4, 5, 6]),

}
### When STH retires, remove indel_file_types!!!
indel_file_types = {
#files
    'pindel.variant_function': ({0: 'var_type_1',
                                 1: 'Gene',
                                 14: 'Sequence',
                                 15: 'Read_Headers',
                                 16: 'Reads',
                                 17: 'Reads2'},
                                [2, 3, 4, 5, 6]),
    'pindel.exonic_variant_function': ({1: 'var_type_2', 2: 'Transcripts'}, [3, 4, 5, 6, 7]),
}

pindel_file_types = {
#files
    'pindel.variant_function': ({0: 'var_type_1',
                                 1: 'Gene',
                                 14: 'Sequence',
                                 15: 'Read_Headers',
                                 16: 'Reads'},
                                [2, 3, 4, 5, 6]),
    'pindel.exonic_variant_function': ({1: 'var_type_2', 2: 'Transcripts'}, [3, 4, 5, 6, 7]),
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
    elif ['GT', 'AD'] == info.keys():
        reads = info['AD'].split(',')
        return reads[0], reads[1], ''
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
        #BCOR,BCOR(NM_001123383:exon8:c.3503-2A>T,NM_001123384:exon7:c.3449-2A>T,NM_017745:exon8:c.3503-2A>T)
        #PHF6(NM_001015877:exon10:c.969-9T>C,NM_032458:exon10:c.969-9T>C),PTEN
        gene_parts=re.split('[(,)]',Gene)
        Gene = gene_parts[0]
        for part in gene_parts:
            if 'NM' in part:
                if Transcripts:
                    Transcripts = ','.join([Transcripts, part])
                else:
                    Transcripts = part
            else:
                if part not in Gene:
                    Gene = ','.join(filter(None, [Gene, part]))
    return Gene, Transcripts

def munge_transcript(data, RefSeqs):
    """
    Return HGVS correct transcript annotations
    Filtered with a preferred transcript list
    NM_006772.1:c.1713G>A
    """
    CODING, PROTEIN = [], []
    transcripts = data.get('Transcripts')
    if transcripts is not None:
        # Split incoming trans, strip the trailing )
        data1 = multi_split(transcripts.replace('(', ':'), ',(')
        #Remove duplicate transcription entries
        data = list(set(data1))
        for d in data:
            codon, prot, protein, coding, txpt = ' ', ' ', ' ', ' ', None
            # Split the actual transcript info which is colon separated
            x = d.split(':')
            #5: ['PRSS1', 'NM_002769', 'exon4', 'c.567T>C', 'p.L189L']
            if len(x)==5:
                gene, txpt, exon, codon, prot = x
            elif len(x)==4:
            #4: ['POLE', 'NM_006231', 'exon25', 'c.2865-4T>-']
                gene, txpt, exon, codon = x
            elif len(x)==3:
            #3: ['RAD50', 'NM_005732', 'c.-38G>A']
                if 'NM' in x[1]:
                    gene, txpt, codon = x
            #3: ['NM_005590','exon5','c.315-4T>-']
                elif 'NM' in x[0]:
                    txpt, exon, codon = x
            elif len(x)==2:
            #2: ['NM_001290310', 'c.*513_*514insATC']
                txpt, codon = x
            elif len(x)==1:
                continue
            else:
                sys.exit("don't know how to parse %s" % d)
            pref_trans = RefSeqs.get(txpt)
            #Want to return None for all values if not pref_trans
            if not pref_trans:
                continue
            code = pref_trans + ':' + codon
            CODING.append(code)
            PROTEIN.append(prot)
    return ' '.join(CODING), ' '.join(PROTEIN)

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
    try:
        if int(data['Var_Reads']) <= 0 and int(data['Ref_Reads']) <= 0:
            freq = 'NA'
        else:
            total = int(data['Var_Reads']) + int(data['Ref_Reads'])
            freq = int(data['Var_Reads']) / float(total)
            freq = "{0:.2f}".format(freq)
    except KeyError:
        freq = 'NA'
    return freq

def munge_ljb_scores(data):
    """
    Parse sift, polyphen and gerp from ljb_all file
    """

    try:
        info = data.get('ljb_Scores').split(',')
        sift = info[0]
        polyphen = info[2]
        mutation_taster = info[8]
        gerp = info[28]

    except AttributeError:
        return -1, -1, -1, -1 

    return polyphen, sift, mutation_taster, gerp

def munge_clinvar(data):
    """
    Combine the clinvar ID into the ClinSig column for linking out in excel spreadsheet
    """
    if data.has_key('ClinVar'):
        CLNALLELEID,CLNDN,CLNDISDB,CLNREVSTAT,CLNSIG=data.get('ClinVar').split(',')
        return ';'.join([CLNALLELEID,CLNSIG]), CLNREVSTAT

    else:
        return -1, -1
    
def munge_cosmic(data):
    if data.has_key('Cosmic'):
        return data.get('Cosmic').split(';')
    else:
        return -1, -1

def largest_variant_reads(output,data):
    """
    return the read info that has the highest variant read
    """
    #Grab the higest read count from pindel
    if data.has_key('Reads2'):
        if data['Reads2'].split(',')[:-1]> data['Reads'].split(',')[:-1]:
            data['Reads']=data['Reads2']

    #Keep the highest variant readcount 
    if output.has_key('Var_Reads') and data.get('Reads'):
        old_ref_read, old_var_read, old_phred = output['Ref_Reads'], output['Var_Reads'], output['Variant_Phred']
        new_ref_read, new_var_read, new_phred = get_reads(data.get('Read_Headers'),data.get('Reads'))
        if int(old_var_read) >  int(new_var_read):
            output['Ref_Reads'], output['Var_Reads'], output['Variant_Phred'] = old_ref_read, old_var_read, old_phred
        else:
            output['Ref_Reads'], output['Var_Reads'], output['Variant_Phred'] = new_ref_read, new_var_read, new_phred
    elif "Var_Reads" not in output.keys() and data.get('Reads'):
        output['Ref_Reads'], output['Var_Reads'], output['Variant_Phred'] = get_reads(data.get('Read_Headers'),data.get('Reads'))

    return output['Ref_Reads'], output['Var_Reads'], output['Variant_Phred']

        
def merge_data(output,data):
    """
    Merge data in certain columns 
    """
    multi_trans_keys=['Transcripts','Gene','var_type_1','var_type_2']
    data_keys=data.keys()
    #If position already in output{}, update certain fields
    for key in multi_trans_keys:
        #if key is already in output dict and in new file
        if key in output.keys() and data.get(key):
            #Make sure the data isn't a duplicate
            if data.get(key) not in output[key]:
                #Make a comma delimited list of info
                output[key]=','.join([output[key],data.get(key)])
            #If position has been seen but this file has no data for this key, set to original data
        elif key in output.keys() and not data.get(key):
            output[key]=output[key]
            #If position has been seen but there was not data for this data, add this file's data 
        elif key not in output.keys() and data.get(key):
            output[key]=data.get(key)
    if output.has_key('Reads') or data.has_key('Reads'):
        output['Ref_Reads'], output['Var_Reads'], output['Variant_Phred'] = largest_variant_reads(output,data)

    multi_trans_keys.extend(['Reads', 'Read_Headers'])
    for k in data_keys:
        if k not in multi_trans_keys:
            output[k]=data.get(k)
    return output

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
        'UW_DEC_p',
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
        'Cosmic_ID',
        'Cosmic_Info',
        'CADD',
        'CLNSIG',
        'CLNREVSTAT',
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

    if args.type == 'PINDEL':
        headers = ['Position'] + variant_headers[3:5] + [
            'Clinically_Flagged',
            'Variant_Type',
            'UW_Freq',
            'Gene',
            'p.',
            'c.',
            'Faves_Y/N',
            'Ref_Reads',
            'Var_Reads',
            'Allele_Frac',
            'Transcripts',
            'UW_Count',
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
        if args.type == 'SNP':
            file_types = snp_file_types
        elif args.type == 'INDEL':
            file_types = indel_file_types
        elif args.type == 'PINDEL':
            file_types = pindel_file_types

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
                continue
                log.warning('no match: %s' % fname)
            if args.strict:
                sys.exit(1)
            continue
        for var_key, data in map_headers(fname, header_ids, var_key_ids):
            if var_key in output:
                output[var_key]=merge_data(output[var_key],data)
            else:
                output[var_key].update(data)
                if output[var_key].has_key('Reads') and not output[var_key].has_key('Var_Reads'):
                    output[var_key]['Ref_Reads'], output[var_key]['Var_Reads'], output[var_key]['Variant_Phred'] = get_reads(data.get('Read_Headers'),data.get('Reads'))

    sort_key = lambda row: [(row[k]) for k in ['chr', 'start', 'stop', 'Ref_Base', 'Var_Base']]
    # # write each row (with all data aggregated), modifying fields as necessary
    for data in sorted(output.values(), key=sort_key):
        variants=[data.get('var_type_2'),data.get('var_type_1')]
        data['Variant_Type'] = ','.join(filter(None, variants))
        data['Gene'], data['Transcripts'] = munge_gene_and_Transcripts(data, RefSeqs)
        data['c.'], data['p.'] = munge_transcript(data, RefSeqs)
        data['Polyphen'], data['Sift'],data['Mutation_Taster'],data['Gerp'] = munge_ljb_scores(data)
        data['Allele_Frac'] = get_allele_freq(data)
        data['CLNSIG'], data['CLNREVSTAT'] = munge_clinvar(data)
        data['Cosmic_ID'],data['Cosmic_Info'] = munge_cosmic(data)
        data['dbSNP_ID'] = data.get('rsid_1') or data.get('rsid_2')
        data['1000g_ALL'] = data.get('1000g_ALL') or -1
        data['1000g_AMR'] = data.get('1000g_AMR') or -1
        data['1000g_SAS'] = data.get('1000g_SAS') or -1
        data['1000g_EAS'] = data.get('1000g_EAS') or -1
        data['1000g_AFR'] = data.get('1000g_AFR') or -1
        data['1000g_EUR'] = data.get('1000g_EUR') or -1
        data['UW_DEC_p'] = data.get('UW_DEC_p') or -1
        data['EXAC'] = data.get('EXAC').split(',')[0] if data.get('EXAC') else -1      
        data['EVS_esp6500_ALL'] = data.get('EVS_esp6500_ALL').split(',')[0] if data.get('EVS_esp6500_ALL') else -1
        data['EVS_esp6500_AA'] = data.get('EVS_esp6500_AA').split(',')[0] if data.get('EVS_esp6500_AA') else -1
        data['EVS_esp6500_EU'] = data.get('EVS_esp6500_EU').split(',')[0] if data.get('EVS_esp6500_EU') else -1
        #CADD is raw score, phred score. We only care about phred
        _, data['CADD'] = split_string_in_two(data.get('CADD'))
        data['ADA_Alter_Splice'],data['RF_Alter_Splice'] = split_string_in_two(data.get('splicing'))
        data['UW_Freq'], data['UW_Count'] = split_string_in_two(data.get('UW_Freq_list'))
        writer.writerow(data)
