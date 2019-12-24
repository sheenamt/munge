"""
Munges the output of different CNV callers into a common, annotated .tsv output
that is usable by the plot_cnv subcommand.
"""

import argparse
import munging.annotation as ann
import numpy as np
import os
import pandas as pd


def build_parser(parser):
    parser.add_argument('cnv_data',
                        help='Path to the raw output from the CNV caller')
    parser.add_argument('package', choices=['contra', 'cnvkit'],
                        help='Software package used to create the CNV_data')
    parser.add_argument('-r', '--refgene',
                        help='Path to the transcript-filtered UCSC RefGene table file (used for annotation)')
    parser.add_argument('-b', '--conifer_baseline', 
                        help='Path to the assay-specific CoNIFER baseline')
    parser.add_argument('-n', '--components_removed', type=int, default=50, 
                        help='Number of components to remove from diagonal matrix during CoNIFER (default: %(default)s)')
    parser.add_argument('-o', '--outfile', 
                        help='Path to the out file (default <prefix>.<package>.CNV_plottable.tsv)')

def parse_contra_file(file_name):
    """Converts a Contra CNATable file to a standard pandas DataFrame for use in plotting."""
    # read in raw output
    df_raw = pd.read_csv(file_name, sep='\t', header=0, dtype=str)
    # create new DataFrame and add relevant columns
    df_plot = pd.DataFrame()
    df_plot['log2'] = df_raw['Adjusted.Mean.of.LogRatio'].astype(float)
    # convert chromosome to standard format ['1', '2', ..., 'X', 'Y', 'GT]
    df_plot['chr'] = df_raw['Chr'].apply(lambda x: ann.chromosomes[x])
    df_plot['start_pos'] = df_raw['OriStCoordinate'].astype(int)
    df_plot['end_pos'] = df_raw['OriEndCoordinate'].astype(int)
    return df_plot

def parse_cnvkit_file(file_name):
    """Converts a CNVkit .cnr file to a standard pandas DataFrame for use in plotting."""
    # read in raw output
    df_raw = pd.read_csv(file_name, sep='\t', header=0, dtype=str)
    # Don't add the antitarget intervals to the plottable df
    df_targets = df_raw[df_raw['gene'] != 'Antitarget']
    # create new DataFrame and add relevant columns
    df_plot = pd.DataFrame()
    df_plot['log2'] = df_targets['log2'].astype(float)
    # convert chromosome to standard format ['1', '2', ..., 'X', 'Y', 'GT]
    df_plot['chr'] = df_targets['chromosome'].apply(lambda x: ann.chromosomes[x])
    df_plot['start_pos'] = df_targets['start'].astype(int)
    df_plot['end_pos'] = df_targets['end'].astype(int)
    return df_plot

def add_annotations(df, genome_tree):
    """
    Adds gene, transcript, and exon number annotations to the standard, plottable
    DataFrame df using the data contained within GenomeIntervalTree genome_tree.

    NOTE: If more than one transcript is found in the genome_tree, annotates with only
    the information from the first transcript found. If row in df covers more than
    one exon in that transcript, labels with only the lowest exon number from that set.
    """
    # note, for performance reasons this should be changed to .apply or something non-iterative
    for i, row in df.iterrows():
        chrom = row['chr']
        start = int(row['start_pos'])
        # searching an interval tree is not inclusive of the endpoint, so increment by one
        end = int(row['end_pos']) + 1
        transcript_list = genome_tree[chrom][start:end]
        if transcript_list:
            # use only the first (and hopefully only) transcript for adding annotations
            t = transcript_list.pop()[2]
            df.at[i, 'gene'] = t.gene
            df.at[i, 'transcript'] = t.id
            # label exons, including those that fall within a UTR
            exons = t.get_exons(start, end, report_utr=False)
            if exons:
                # use only the first (and hopefully only) exon number for each interval
                df.at[i, 'exon'] = str(exons[0].number)
        # if no hits from the genome tree, annotate as intergenic
        else:
            df.at[i, 'gene'] = 'intergenic'

def run_conifer(sample_df, baseline_df, components_removed):
    return sample_df

def action(args):
    # import and parse the data
    if args.package == 'contra':
        df = parse_contra_file(args.cnv_data)
    elif args.package == 'cnvkit':
        df = parse_cnvkit_file(args.cnv_data)
    else:
        # should never hit here
        raise ValueError("Improper package specified as argument")

    # add annotations to the parsed data
    gt = ann.GenomeIntervalTree.from_table(args.refgene)
    add_annotations(df, gt)

    # apply CoNIFER if requested
    if args.conifer_baseline:
        # read in the conifer baseline
        baseline_df = pd.read_feather(args.conifer_baseline)
        baseline_df = baseline_df.set_index(['chr', 'start_pos', 'end_pos'])
        # drop duplicate probes from df
        df = df.set_index(['chr', 'start_pos', 'end_pos'])
        df = df.loc[~df.index.duplicated(keep='first')]
        # merge sample into baseline using index
        baseline_df = pd.merge(baseline_df, df[['log2']], how='left', left_index=True, right_index=True)
        # impute missing log2 values in sample column with probe medians
        baseline_df['log2'] =  baseline_df['log2'].fillna(baseline_df['probe_median'])
        # transform
        log2_matrix = baseline_df[baseline_df.columns[2:]].values
        U, S, Vt = np.linalg.svd(log2_matrix,full_matrices=False)
        new_S = np.diag(np.hstack([np.zeros([args.components_removed]),S[args.components_removed:]]))
        transformed_log2s = np.dot(U, np.dot(new_S, Vt))
        # merge the transformed log2 column for the sample into the plottable df
        baseline_df['conifer_log2'] = transformed_log2s[:,-1]
        df = df.merge(baseline_df[['conifer_log2']], how='inner', left_index=True, right_index=True)
        out_columns=['chr', 'start_pos', 'end_pos', 'log2', 'conifer_log2', 'gene', 'transcript', 'exon']
    else:
        out_columns=['chr', 'start_pos', 'end_pos', 'log2', 'gene', 'transcript', 'exon']

    # set the outfile path
    if args.outfile:
        out_path = args.outfile
    else:
        dir_name = os.path.dirname(args.cnv_data)
        base_name = os.path.basename(args.cnv_data)
        prefix = base_name.split('.')[0]
        out_file = "{}.{}.CNV_plottable.tsv".format(prefix, args.package)
        out_path = os.path.join(dir_name, out_file)
    
    # make chr column sortable with natural sorting order
    df = df.reset_index()
    df['chr'] = pd.Categorical(df['chr'], categories=ann.chromosome_sort_order, ordered=True)
    # sort by chromosome, then position
    print(df.dtypes)
    df = df.sort_values(['chr', 'start_pos']) 
    # save the data
    df.to_csv(out_path, sep='\t', index=False, columns=out_columns)