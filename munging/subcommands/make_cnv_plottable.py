import argparse
import os
import pandas as pd

def build_parser(parser):
    parser.add_argument('cnv_data',
                        help='Path to the raw output from the CNV caller')
    parser.add_argument('package', choices=['contra', 'cnvkit'],
                        help='Software package used to create the CNV_data')
    parser.add_argument('-r', '--refgene',
                        help='Path to the transcript-filtered USCS RefSeq table file')
    parser.add_argument('-o', '--outfile', 
                        help='Path to the out file (default <prefix>.<package>.CNV_plottable.tsv)')

def parse_contra_file(file_name, refseq_table):
    df_raw = pd.read_csv(file_name, sep=',', header=0, dtype=str)
    df_plot = pd.DataFrame()
    df_plot['log2'] = df_raw['Adjusted.Mean.of.LogRatio']
    df_plot['chr'] = df_raw['Chr']
    df_plot['start_pos'] = df_raw['OriStCoordinate']
    df_plot['end_pos'] = df_raw['OriEndCoordinate']
    gene_symbols = df_raw['Gene.Sym'].str.split(':', n = 1, expand = True)
    df_plot['gene'] = gene_symbols[0]       # in future, change this to using munge annotation
    df_plot['transcript'] = gene_symbols[1] # in future, change this to using munge annotation
    df_plot['exon'] = df_raw['Exon.Number']    # in future, change this to using munge annotation

    return df_plot

def parse_cnkvit_file(file_name, refseq_table):
    # TO DO
    df_raw = pd.read_csv(file_name, sep='/t', header=0)
    df_plot = pd.DataFrame()
    df_plot['log2'] = [0]
    df_plot['chr'] = [0]
    df_plot['start_pos'] = [0]
    df_plot['end_pos'] = [0]
    df_plot['gene'] = [0]
    df_plot['transcript'] = [0]
    df_plot['exon'] = [0]

    return df_plot

def action(args):
    # import and parse the data
    if args.package == 'contra':
        df = parse_contra_file(args.data, args.refgene)
    elif args.package == 'cnvkit':
        df = parse_cnvkit_file(args.data, args.refgene)
    else:
        # should never hit here
        raise ValueError("Improper type specified as argument")

    # set the out file path
    if args.outfile:
        out_path = args.file
    else:
        dir_name = os.path.dirname(args.cnv_data)
        base_name = os.path.basename(args.cnv_data)
        prefix = base_name.split('.')[0]
        out_file = "{}.{}.CNV_plottable.tsv".format(prefix, args.package)
        out_path = os.path.join(dir_name, base_name)
    # save the data
    df.to_csv(out_path, sep='\t', index=False)