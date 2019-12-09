"""
Munges the output of different CNV callers into a common, annotated .tsv output
that is usable by the plot_cnv subcommand.
"""

import argparse
import os
import pandas as pd
import munging.annotation as ann

def build_parser(parser):
    parser.add_argument('cnv_data',
                        help='Path to the raw output from the CNV caller')
    parser.add_argument('package', choices=['contra', 'cnvkit'],
                        help='Software package used to create the CNV_data')
    parser.add_argument('-r', '--refgene',
                        help='Path to the transcript-filtered UCSC RefGene table file (used for annotation)')
    parser.add_argument('-o', '--outfile', 
                        help='Path to the out file (default <prefix>.<package>.CNV_plottable.tsv)')

def parse_contra_file(file_name):
    """Converts a Contra CNATable file to a standard pandas DataFrame for use in plotting."""
    df_raw = pd.read_csv(file_name, sep='\t', header=0, dtype=str)
    df_plot = pd.DataFrame()
    df_plot['log2'] = df_raw['Adjusted.Mean.of.LogRatio']
    # convert chromosome to standard format ['1', '2', ..., 'X', 'Y', 'GT]
    df_plot['chr'] = df_raw['Chr'].apply(lambda x: ann.chromosomes[x])
    # make chr column sortable with natural sorting order
    df_plot['chr'] = pd.Categorical(df_plot['chr'], categories=ann.chromosome_sort_order, ordered=True)
    df_plot['start_pos'] = df_raw['OriStCoordinate']
    df_plot['end_pos'] = df_raw['OriEndCoordinate']
    return df_plot

def parse_cnvkit_file(file_name):
    """Converts a CNVkit .cnr file to a standard pandas DataFrame for use in plotting."""
    df_raw = pd.read_csv(file_name, sep='\t', header=0, dtype=str)
    # Don't add the antitarget intervals to the plottable df
    df_targets = df_raw[df_raw['gene'] != 'Antitarget']
    df_plot = pd.DataFrame()
    df_plot['log2'] = df_targets['log2']
    # convert chromosome to standard format ['1', '2', ..., 'X', 'Y', 'GT]
    df_plot['chr'] = df_targets['chromosome'].apply(lambda x: ann.chromosomes[x])
    # make chr column sortable with natural sorting order
    df_plot['chr'] = pd.Categorical(df_plot['chr'], categories=ann.chromosome_sort_order, ordered=True)
    df_plot['start_pos'] = df_targets['start']
    df_plot['end_pos'] = df_targets['end']
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

    # sort by chromosome, then position
    df = df.sort_values(['chr', 'start_pos']) 

    # set the outfile path
    if args.outfile:
        out_path = args.outfile
    else:
        dir_name = os.path.dirname(args.cnv_data)
        base_name = os.path.basename(args.cnv_data)
        prefix = base_name.split('.')[0]
        out_file = "{}.{}.CNV_plottable.tsv".format(prefix, args.package)
        out_path = os.path.join(dir_name, out_file)
    # save the data
    df.to_csv(out_path, sep='\t', index=False)