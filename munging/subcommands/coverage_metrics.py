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
        df = pd.read_csv(f, comment='#', delimiter='\t',header=None,usecols=[0,1,4,5],names=headers)
        df['POS']=df['START']+(df['1_BASED_POSITION']-1)

        #remove all rows that pass the thresholds
        failing_bases = df[df['READS']<=args.threshold]

        #Now, groupy by START, mean READS, and label interval with first and last POS
        data=failing_bases.groupby(['CHROM','START'],as_index=False).agg({'READS':['mean','median'],
                                                                          'POS':['min','max']})
        
        #SEt interval to FALSE/TRUE if SNP or NOT 
        data.loc[(data['POS']['min'] == data['POS']['max']),'INTERVAL']=False
        data.loc[(data['POS']['min'] < data['POS']['max']),'INTERVAL']=True

        #NOW, ANNOTATE:
        
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
    print(output)
    sys.exit()
    out_fieldnames=['Position','Gene','Gene_Region','Transcripts','Mean Coverage']
    writer = csv.DictWriter(args.outfile, extrasaction='ignore',fieldnames=out_fieldnames, delimiter='\t')
    writer.writeheader()
    writer.writerows(output) 


        data.loc[
    #now that we hace the read depth info, and the region, we need to annotate it
    #first, only process real chromosoms

    #if min==max, we're dealing with a SNP and therfore only 1 postions
#    ,data['POS']['min'])    
    data=data.apply(define_interval)
    sys.exit()
    #otherwise, we're dealing wiht an interval

    #rows = data.T.to_dict().values()
    # #annotate failing intervals

    for row in rows:
        print(row)
        sys.exit()
        try:
            if chrm:
                chr1 = 'chr'+str(chromosomes[row['CHROM']])
            else:
                chr1 = str(chromosomes[row['CHROM']])
        except KeyError:
            continue

        sys.exit()
        #Setup the variables to be returned
        gene1, region, transcripts=[],[],[]
        # each segment is assigned to a gene if either the
          # start or end coordinate falls within the feature boundaries.
        snp=False
#        if int(row['START']) == int(row['END']) :
        if int(row['POS']) == int(row['END']) :
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
    print(output)
    sys.exit()
    out_fieldnames=['Position','Gene','Gene_Region','Transcripts','Mean Coverage']
    writer = csv.DictWriter(args.outfile, extrasaction='ignore',fieldnames=out_fieldnames, delimiter='\t')
    writer.writeheader()
    writer.writerows(output) 


    
#         # # First row has start that is preceded by different start
#         # fst=df.index[df['START'] and not df['START'].shift(1).fillna(False)]

#         # print(fst)

#         print((df.loc[df['CHROM']==df['CHROM'] ]) and (df.loc[df['POS']==df['POS'] ])) # and (df['READS']<=args.threshold))

# #        varscan_format_variants.loc[(varscan_format_variants['chrom'] == chrom) & (varscan_format_variants['varscan_start'] == pos_start), 'Valid_Reads'] = depth
#         print(fst)


#         # first row is a True preceded by a False
#         #" & ~ "== "and not"
#         fst = df.index[df['tag'] & ~ df['tag'].shift(1).fillna(False)]
#         #fst = df.index[df['POS'] &  df['POS'].shift(1).fillna(False)]
#         # last row is a True followed by a False
#         lst = df.index[df['tag'] & ~ df['tag'].shift(-1).fillna(False)]

#         print('fst:', fst)
#         print('lst:',lst)
#         # create list of first and last in interval
#         pr = [(i, j) for i, j in zip(fst, lst)] # if j > i + 4]
#         for i,j in zip(fst,lst):
#             print(i)
#             print(j)
#         sys.exit()
#         print('pr:', pr)
#         #Empty dict for adding data that will need annotations
#         #parse original datafram into intervals 
#         data=pd.DataFrame(columns=['CHROM','START','END','Mean Coverage'])
#         for i, j in pr:
#             new_df=df.loc[i:j]
#             mean_data = int(new_df['READS'].mean())
#             chrom=new_df.iloc[0]['CHROM']
#             start=new_df.iloc[0]['POS']
#             end=new_df.iloc[-1]['POS']
#             data=data.append({'CHROM':str(chrom),'START':str(start),'END':str(end),'Mean Coverage':mean_data}, ignore_index=True)
 
