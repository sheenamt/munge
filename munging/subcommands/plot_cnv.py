"""
ADD DOCUMENTATION
"""

import argparse
import os
import pandas as pd
# needed to avoid Xwindows backend when using python 2.7 (which will cause an error)
import matplotlib
matplotlib.use('Agg')   # must call before importing pyplot with python 2.7
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from natsort import natsorted


def build_parser(parser):
    parser.add_argument('cnv_data',
                        help='Path to the CNV_plottable.tsv file')
    parser.add_argument('outfile', 
                        help='Path to out file')
    parser.add_argument('-t', '--min_log_ratio', type=float, default=0.5, 
                        help='Minimum abs(log ratio) for printing gene names (default: %(default)s)')
    parser.add_argument('-w', '--window_size', type=int, default=20, 
                        help='Window size for rolling median (default: %(default)s)')
    parser.add_argument('--title', type=str, 
                        help='Title (default: infer from cnv_data file name)')

def load_cnv_data(file_path, window_size):
    column_data_types = {'log2' : 'float64',
                         'chr' : 'str',
                         'start_pos' : 'int64',
                         'end_pos' : 'int64',
                         'gene' : 'str',
                         'transcript' : 'str',
                         'exon' : 'str'
                         }
    df = pd.read_csv(file_path, sep='\t', header=0, dtype=column_data_types)
    df['mean_pos'] = df[['start_pos', 'end_pos']].mean(axis=1).astype(int)
    df['rolling_median'] = df['log2'].rolling(window_size, min_periods=1).median()
    return df

def extract_plot_title(file_path):
    """Returns the 'genetics string' of the data file"""
    base_name = os.path.basename(file_path)
    return base_name.split('.')[0]

def flag_genes(df, min_log_ratio):
    """
    Flags genes that are outside the bounds specified by min_log_ratio

    Genes are flagged if their most extreme 'rolling_median' is outside the bounds
    defined by [-min_log_ratio, min_log_ratio]

    Returns a dict of <gene>:(<index_of_rolling_median> ,<value_rolling_median>) 
    """
    flagged_genes = {}
    for gene in df['gene'].unique():
        if gene == 'intergenic':
            continue
        df_gene = df[df['gene']==gene]
        # find the max median and its index, and flag gene if condition met
        max_idx = df_gene['rolling_median'].idxmax()
        max_median = df_gene.loc[max_idx]['rolling_median']
        if max_median > min_log_ratio:
            flagged_genes[gene] = (max_idx, max_median)
            continue
        # find the min median and its index, and flag gene if condition met
        min_idx = df_gene['rolling_median'].idxmin()
        min_median = df_gene.loc[min_idx]['rolling_median']
        if min_median < -1 * min_log_ratio:
            flagged_genes[gene] = (min_idx, min_median)
        
    return flagged_genes

def plot_main(pdf, df, title, min_log_ratio):
    """Saves to the pdf a plot of every point across all chromosomes present in df, flagging genes above min_log_ratio"""
    fig, ax =  plt.subplots(figsize=(11, 8.5))
    
    # set plot scale
    Y_SCALE = 2
    plt.ylim((-1 * Y_SCALE, Y_SCALE))

    # plot entries left to right from chrom 1 to chrom X
    chromosomes = natsorted(df['chr'].unique())

    x_tick_values = []
    v_line_coords = []
    for chrom in chromosomes:
        df_chrom = df[df['chr'] == chrom]
        indices = df_chrom.index.values
        ax.scatter(indices, df_chrom['log2'], label=chrom, marker=',', s=0.2)
        x_tick_values.append(np.median(indices))
        v_line_coords.append(np.max(indices))

    # flag genes that are above or below log ratio threshold
    flagged_genes = flag_genes(df, min_log_ratio)

    # add labels to plot for genes that are flagged (shift labels if needed to appear on chart)
    shifted_labels = {}
    for gene, coord in flagged_genes.items():
        label_shifted = False
        if coord[1] > Y_SCALE:
            coord = (coord[0], 0.97 * Y_SCALE)
            label_shifted = True
        elif coord[1] < -1 * Y_SCALE:
            coord = (coord[0], -0.99 * Y_SCALE)
            label_shifted = True        
        a = plt.annotate(gene, coord, fontsize=8)
        if label_shifted:
            shifted_labels[gene] = a

    # set plot ticks and tick labels
    plt.tick_params(right=True, top=True)
    plt.xticks(x_tick_values, chromosomes, fontsize=8)
    plt.yticks(fontsize=10)
    plt.locator_params(axis='y', nbins=4)
    
    # add vertical lines (skip line after last chrom)
    for coord in v_line_coords[0:-1]:
        plt.axvline(x=coord, alpha=0.1, color='black', linewidth=1)

    # add horizontal lines
    plt.grid(axis='y', which='major', linestyle='solid', linewidth=1)
    for coord in [min_log_ratio * -1, min_log_ratio]:
        plt.axhline(y=coord, alpha=0.5, color='black', linewidth=1, linestyle='dotted')

    # label axes
    plt.xlabel('Chromosome', fontsize=12)
    plt.ylabel('Adjusted Mean of Log Ratio', fontsize=12)

    # set title and legend
    plt.title(title)
    plt.legend(fontsize=10, 
               markerfirst=False, 
               markerscale=20, 
               bbox_to_anchor=(1.02, 0.5), 
               borderaxespad=0.2, 
               loc='center left')
    
    # save main figure
    pdf.savefig()

    # if logs2 were cutoff by [-2,2] plot, add second plot covering full range
    min_log = df['log2'].min()
    max_log = df['log2'].max()
    if min_log < -1 * Y_SCALE or max_log > Y_SCALE:
        plt.ylim((min_log - 0.1, max_log + 0.1))
        plt.title(title + ' (Plot 2)')
        #ax.yaxis.set_major_locator(AutoLocator())

        # move labels that were shifted in original plot
        for gene, label in shifted_labels.items():
            label.set_y(flagged_genes[gene][1])
            #label.remove()
            #plt.annotate(gene, flagged_genes[gene], fontsize=8)

        # save second figure
        pdf.savefig()

    plt.close()

def plot_gene(pdf, df_gene):
    """
    Saves to pdf a plot of every log2 in df_gene, colored and labeled by exon
    """
    fig, ax =  plt.subplots(figsize=(11, 8.5))
    
    # set plot scale
    Y_SCALE = 2
    plt.ylim((-1 * Y_SCALE, Y_SCALE))

    # plot non-exonic entries
    df_introns = df_gene[pd.isnull(df_gene['exon'])]
    ax.scatter(df_introns['mean_pos'], df_introns['log2'], marker='o', color='black', alpha=0.7, s=3)
    
    # parameters for y-coordinate of exon label
    exon_labels=[]
    gene_median = df_gene['log2'].median()
    if gene_median > 0:
        # labels go below the median
        offset_factor = -1
    else:
        # labels go above the median
        offset_factor = 1
    max_label_offset = 1
    # stagger the labels into 3 rows to reduce overlapping text
    offset_step_size = 0.1
    max_offset_steps = 3

    # plot exonic entries
    exons = df_gene['exon'].dropna().unique()
    for i, exon in enumerate(exons):
        # plot entries for exon
        df_exon = df_gene[df_gene['exon'] == exon]
        s = ax.scatter(df_exon['mean_pos'], df_exon['log2'], label=exon, marker='o', s=6)
        # add label for exon
        color = s.get_facecolor()[0] 
        label_x = df_exon['mean_pos'].mean()
        label_y = gene_median + offset_factor * (max_label_offset - offset_step_size * (i % max_offset_steps))
        a = plt.annotate(exon, (label_x, label_y), fontsize=10, color=color)
        exon_labels.append(a)

    # set plot ticks and tick labels
    plt.tick_params(right=True, top=True)
    plt.xticks(fontsize=10)
    plt.locator_params(axis='x', nbins=4)
    plt.ticklabel_format(axis='x', style='plain', useOffset=False)
    plt.yticks(fontsize=10)
    plt.locator_params(axis='y', nbins=4)
    
    # add horizontal lines
    plt.grid(axis='y', which='major', linestyle='solid', linewidth=1)

    # label axes
    plt.xlabel('Position', fontsize=12)
    plt.ylabel('Adjusted Mean of Log Ratio', fontsize=12)

    # set title
    title_parts = df_gene.iloc[0][['chr', 'gene', 'transcript']]
    plt.title("Chromosome {}: {}:{}".format(*title_parts))
    
    # save main figure
    pdf.savefig()

    # if log2s were cutoff by [-2,2] plot, add second plot covering full range
    min_log = df_gene['log2'].min()
    max_log = df_gene['log2'].max()
    if min_log < -1 * Y_SCALE or max_log > Y_SCALE:
        plt.ylim((min_log - 0.1, max_log + 0.1))
        plt.title("Chromosome {}: {}:{} (Plot 2)".format(*title_parts))

        # move exon labels to account for new scale
        y_mid = (max_log + min_log) / 2
        if gene_median < y_mid:
            # labels go above median log
            new_factor = (max_log - min_log) / (Y_SCALE * 2)
        else:
            # labels go below median log
            new_factor = -1 * (max_log - min_log) / (Y_SCALE * 2)

        for a in exon_labels:
            old_y = a.xy[1]
            new_y = (old_y - gene_median) * (new_factor / offset_factor) + gene_median
            a.set_y(new_y)
        # save second figure
        pdf.savefig()

    plt.close()

def action(args):
    # load the data
    df = load_cnv_data(args.cnv_data, args.window_size)

    # infer title for main plot if needed
    if not args.title:
        args.title = extract_plot_title(args.cnv_data)

    # create plots and save to pdf-formatted outfile
    with PdfPages(args.outfile) as pdf:

        plot_main(pdf, df, args.title, args.min_log_ratio)

        for gene in df['gene'].unique():
            if gene == 'intergenic':
                continue
            gene_df = df[df['gene']==gene]
            plot_gene(pdf, gene_df)