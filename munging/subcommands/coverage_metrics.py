"""
Parse bedtools output for intervals not meeting given threshold, annotate with RefSeq data

Usage:

munge coverage_metrics /path/to/bedtools/output /path/to/refseq -o /path/to/output
"""
import sys
import csv
import pandas as pd
from munging.utils import Opener
from munging.annotation import chromosomes,GenomeIntervalTree, UCSCTable, define_transcripts

def build_parser(parser):
    parser.add_argument('target_coverage',
                        help='Bedtools output for base level coverage')
    parser.add_argument('refgene', 
                        help='RefGene file')
    parser.add_argument('-t','--threshold', type=int, default=50,
                        help='Threshold for outputing failed coverage')
    parser.add_argument('-o', '--outfile', type=Opener('w'), metavar='FILE',
                        default=sys.stdout, help='output file')

def action(args):
    #read in refgene into interval tree
    #Create interval tree of introns and exons,  grouped by chr
    exons = GenomeIntervalTree.from_table(open(args.refgene, 'r'), parser=UCSCTable.REF_GENE, mode='exons')

    #try to match chr string in reference file and input data
    chrm=False
    if 'chr' in exons.keys()[0]:
        chrm = True

    output = []

    #read in bedfile
    with open(args.target_coverage, 'rU')  as f:
        headers=['CHROM','START','END','GENE','1_BASED_POSITION','READS']
        reader = pd.read_csv(f, comment='#', delimiter='\t',header=None,usecols=[0,1,4,5],names=headers)
        reader['POS']=reader['START']+(reader['1_BASED_POSITION']-1)
        
        #parse into intervals not meeting coverage threshold
        failing_bases=reader[reader['READS']<=args.threshold]
        df = reader
        df['tag'] = df['READS']<=args.threshold

        # first row is a True preceded by a False
        fst = df.index[df['tag'] & ~ df['tag'].shift(1).fillna(False)]
        # last row is a True followed by a False
        lst = df.index[df['tag'] & ~ df['tag'].shift(-1).fillna(False)]

        # create list of first and last in interval
        pr = [(i, j) for i, j in zip(fst, lst)] # if j > i + 4]
        
        #Empty dict for adding data that will need annotations
        #parse original datafram into intervals 
        data=pd.DataFrame(columns=['CHROM','START','END','Mean Coverage'])
        for i, j in pr:
            new_df=df.loc[i:j]
            mean_data = int(new_df['READS'].mean())
            chrom=new_df.iloc[0]['CHROM']
            start=new_df.iloc[0]['POS']
            end=new_df.iloc[-1]['POS']
            data=data.append({'CHROM':str(chrom),'START':str(start),'END':str(end),'Mean Coverage':mean_data}, ignore_index=True)
                    
    data=data.drop_duplicates()
    #Convert to a dictionary for processing clearly
    rows = data.T.to_dict().values()

    # #annotate failing intervals

    for row in rows:
        #only normal chr are process, GL amd MT are ignored
        try:
            if chrm:
                chr1 = 'chr'+str(chromosomes[row['CHROM']])
            else:
                chr1 = str(chromosomes[row['CHROM']])
        except KeyError:
            continue
        #Setup the variables to be returned
        gene1, region, transcripts=[],[],[]
        # each segment is assigned to a gene if either the
        # start or end coordinate falls within the feature boundaries.
        snp=False
        if int(row['START']) == int(row['END']) :
            snp=True
            chrm_exons=exons[chr1].search(int(row['START']))
        else:
            chrm_exons=exons[chr1].search(int(row['START']), int(row['END']))
            chrm_start=exons[chr1].search(int(row['START']))
            chrm_stop=exons[chr1].search(int(row['END']))

        #Usual case: both start and stop are in a coding region
        if chrm_exons:
            gene1, region, transcripts=define_transcripts(chrm_exons)
        # #But if start isn't in coding, but stop is, currently an error
        elif chrm_stop and not chrm_start:
            print 'hit this case'
            sys.exit()
            gene1, region, transcripts=define_transcripts(chrom_stop)
        #Otherwise if neither start nor stop are in coding, label everything as intergenic
        else:
            gene1=['Intergenic',]
            region=['Intergenic',]
            transcripts=[]

        row['Gene'] =';'.join(str(x) for x in set(gene1))
        row['Gene_Region']=';'.join(str(x) for x in set(region))
        row['Transcripts']=';'.join(str(x) for x in set(transcripts))

        if snp:
            row['Position']=str(chr1)+':'+str(row['START'])
        else:
            row['Position']=str(chr1)+':'+str(row['START'])+'-'+str(row['END'])
        output.append(row)

    out_fieldnames=['Position','Gene','Gene_Region','Transcripts','Mean Coverage']
    writer = csv.DictWriter(args.outfile, extrasaction='ignore',fieldnames=out_fieldnames, delimiter='\t')
    writer.writeheader()
    writer.writerows(output) 
