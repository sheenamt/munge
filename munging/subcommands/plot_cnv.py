"""
Creates annotated plots of adjusted log2 ratios from CNV callers and saves the
plots to a PDf.

PLOTS contra and conifer on the same set of axes
"""

import argparse
import os
import sys
import pandas as pd
import numpy as np
# needed to avoid Xwindows backend when using python 2.7 (which will cause an error)
import matplotlib
matplotlib.use('Agg')   # must call before importing pyplot with python 2.7
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
from natsort import natsorted
from munging.annotation import UCSCTable, Transcript

# global parameter for the preferred limits to the y-axis [-2, 2]
y_scale = 2

def build_parser(parser):
    parser.add_argument('cnv_data',
                        help='Path to the CNV_plottable.tsv file')
    parser.add_argument('-o', '--outfile', default='plot_cnv.pdf',
                        help='Path to out file (default: %(default)s)')
    parser.add_argument('-t', '--min_log_ratio', type=float, default=0.5, 
                        help='Minimum abs(log ratio) for printing gene names (default: %(default)s)')
    parser.add_argument('-w', '--window_size', type=int, default=20, 
                        help='Window size for rolling median (default: %(default)s)')
    parser.add_argument('-r', '--refgene',
                        help='Path to the transcript-filtered USCS RefSeq table file (used to add IGV-like figures to gene-level plots)')
    parser.add_argument('--title', type=str, 
                        help='Title (default: infer from cnv_data file name)')

def load_cnv_data(file_path):
    """Returns a DataFrame loaded from file path, adding columns for mean position and rolling median based on window size"""
    column_data_types = {'log2' : 'float64',
                         'conifer' : 'float64',
                         'chr' : 'str',
                         'start_pos' : 'int64',
                         'end_pos' : 'int64',
                         'gene' : 'str',
                         'transcript' : 'str',
                         'exon' : 'str'
                         }
    df = pd.read_csv(file_path, sep='\t', header=0, dtype=column_data_types)
    df['mean_pos'] = df[['start_pos', 'end_pos']].mean(axis=1).astype(int)
    return df

def filled_rolling_median(series, window_size):
    """
    Returns the median of a center-aligned rolling window of size `window_size` for `series`.
    
    Rather than returning NaN at the beginning and end of the series, uses a left-aligned or right-aligned window.
    If the series is shorter than the window size, shrinks window to series length.
    """
    # shrink the window if it's larger than the length of the series
    if len(series) < window_size:
        window_size = len(series)

    # fill nan values (will only apply to probes not covered by conifer)
    no_nan = series.fillna(0)

    # create left, centered, and right rolling medians
    right = no_nan.rolling(window_size, center=False).median()
    center = no_nan.rolling(window_size, center=True).median()
    left = no_nan.rolling(window_size, center=False).median().shift(-1 * window_size)
    
    # fill the gaps in center with values from either end
    filled = center.fillna(right).fillna(left)

    return filled

def extract_plot_title(file_path):
    """Returns the 'genetics string' of the data file"""
    base_name = os.path.basename(file_path)
    return base_name.split('.')[0]

def flag_genes(df, min_log_ratio, rolling_window_size):
    """
    Flags genes that are outside the bounds specified by min_log_ratio

    Genes are flagged if their most extreme 'rolling_median' is outside the bounds defined by 
    -min_log_ratio, min_log_ratio] or if four consecutive conifer points are outside said bounds and they deviate
    strongly from the gene median.

    Returns a dict of <gene>:(<index_of_gene_extrema>, <value_gene_extrema) 
    """
    flagged_genes = {}
    data_column_name = df.columns[0]  # allows for flagging conifer or non-conifer log2s

    for gene in df['gene'].unique():
        # don't flag intergenic entries
        if gene == 'intergenic':
            continue
        # create a subslice of entries for that gene
        df_gene = df[df['gene']==gene].copy()

        # find extreme points of the gene
        max_x = df_gene[data_column_name].idxmax()
        max_y = df_gene.loc[max_x][data_column_name]
        min_x = df_gene[data_column_name].idxmin()
        min_y = df_gene.loc[min_x][data_column_name]

        # calculate a rolling median at each probe
        df_gene['rolling_median'] = filled_rolling_median(df_gene[data_column_name], rolling_window_size)

        # find the max median and its index, and flag gene if condition met
        max_median = df_gene['rolling_median'].max()
        if max_median > min_log_ratio:
            flagged_genes[gene] = (max_x, max_y)
            continue

        # find the min median and its index, and flag gene if condition met
        min_median = df_gene['rolling_median'].min()
        if min_median < -1 * min_log_ratio:
            flagged_genes[gene] = (min_x, min_y)
            continue

        # if conifer, flag if at least 4 consecutive probes are above or below threshold AND
        # the distance between the gene median and the probes is larger than the threshold
        if data_column_name == 'conifer':
            CONIFER_WINDOW = 4
            median = df_gene['conifer'].median()
            # flag probes
            df_gene['flagged'] = [True if min(abs(x), abs(x - median)) > min_log_ratio else False for x in df_gene['conifer']]
            # count flagged probes in rolling window of 4
            df_gene['rolling_flagged'] = df_gene['flagged'].rolling(window=CONIFER_WINDOW).sum()
            # flag gene if condition met
            if df_gene['rolling_flagged'].max() >= CONIFER_WINDOW:
                label_x = df_gene['rolling_flagged'].idxmax()
                label_y = df_gene.loc[label_x]['conifer']
                flagged_genes[gene] = (label_x, label_y)
                continue
        
    return flagged_genes

def create_transcript_subplot(axis, transcript):
    """Creates an IGV-like representation of transcript on axis"""
    left, right = axis.get_xlim()

    # add horizontal line from transcription start to end
    axis.hlines(0.5, left, right, color='b', linewidth=1)

    # add arrows to indicate direction of transcription
    if transcript.strand == '+':
        arrow_text = '>'
    else:
        arrow_text = '<'
    
    NUM_ARROWS = 15     # hard coded parameter for number of arrows to add
    arrow_step = (right - left)*1.0 / NUM_ARROWS

    for i in range(NUM_ARROWS):
        arrow_x = left + (i + 0.5) * arrow_step
        axis.text(arrow_x, 0.48, arrow_text, color='b', fontsize=10, ha='right', va='center')

    # add boxes for exons
    exons = transcript.get_exons(report_utr=False)
    boxes = []

    for e in exons:
        # if the entire exon is non-coding, make a small box for whole exon and continue
        if e.cd_start is None or e.cd_end is None:
            length = e.end - e.start + 1
            r = patches.Rectangle((e.start, 0.375), width=length, height=0.25, color='b')
            boxes.append(r)
            continue
        # if the start of the exon is non-coding, make a small box for non-coding region at start
        if e.cd_start > e.start:
            length = e.cd_start - e.start + 1
            r = patches.Rectangle((e.start, 0.375), width=length, height=0.25, color='b')
            boxes.append(r)
        # if the end of the exon is non-coding, make a small box for non-coding region at end
        if e.cd_end < e.end:  
            length = e.end - e.cd_end + 1
            r = patches.Rectangle((e.cd_end, 0.375), width=length, height=0.25, color='b')
            boxes.append(r)

        # make a large box for coding region of exon
        length = e.cd_end - e.cd_start + 1
        r = patches.Rectangle((e.cd_start, 0.25), width=length, height=0.5, color='b')
        boxes.append(r)

    for b in boxes:
        axis.add_patch(b)

def create_gene_subplot(axis, df_subplot, data_column):
    """
    Plots a gene-level subplot to axis using the information contained within df_subplot and data_column
    
    Returns a list of all exon labels created by the plot (for use in the event that the plot needs to be rescaled)
    """
    df_introns = df_subplot[pd.isnull(df_subplot['exon'])]
    covered_exons = df_subplot['exon'].dropna().unique()

    # plot non-exonic entries in black
    axis.scatter(df_introns['mean_pos'], df_introns[data_column], marker='o', color='black', alpha=0.7, s=3)

    # parameters for y-coordinate of exon label
    exon_labels=[]
    gene_median = df_subplot[data_column].median()
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

    # plot exonic entries, each in their own color
    exon_positions = [] # needed for box and whiskers plot
    exon_vectors = []   # needed for box and whiskers plot

    for i, exon in enumerate(covered_exons):
        # plot entries for exon
        df_exon = df_subplot[df_subplot['exon'] == exon]
        s = axis.scatter(df_exon['mean_pos'], df_exon[data_column], label=exon, marker='o', s=12)

        # add label for exon
        color = s.get_facecolor()[0] 
        label_x = df_exon['mean_pos'].mean()
        label_y = gene_median + offset_factor * (max_label_offset - offset_step_size * (i % max_offset_steps))
        a = axis.text(label_x, label_y, exon, fontsize=10, color=color, ha='center')
        exon_labels.append(a)
        # add exon data for box and whiskers plot
        exon_positions.append(df_exon['mean_pos'].median())
        exon_vectors.append(df_exon[data_column].values)

    # add a box and whiskers around each exon
    min_x = df_subplot['mean_pos'].min()
    max_x = df_subplot['mean_pos'].max()
    box_width = (max_x - min_x) / 80
    median_props = dict(linestyle='-', linewidth=3, color='black')
    axis.boxplot(exon_vectors,
            positions=exon_positions,
            manage_xticks=False,
            widths=box_width,
            autorange=True,
            medianprops= median_props,
            showfliers=False)

    # set plot limits, ticks, and tick labels
    x_pad = int((max_x - min_x) * 0.04)
    if x_pad > 0:
        axis.set_xlim(min_x - x_pad, max_x + x_pad)
    axis.set_ylim((-1 * y_scale, y_scale))
    axis.locator_params(axis='both', nbins=4)
    axis.tick_params(right=True, top=True, labelsize=10)
    axis.ticklabel_format(axis='x', style='plain', useOffset=False)
    
    # add horizontal lines
    axis.grid(axis='y', which='major', linestyle='solid', linewidth=1)

    return exon_labels

def plot_gene(pdf, df_gene, transcript=None):
    """
    Saves to pdf a plot of every log2 in df_gene, colored and labeled by exon

    If a Transcript object is provided for transcript, also creates an IGV-like representation of
    the gene at the bottom of the plot.

    If df_gene contains a conifer column, adds a subplot to the conifer data
    """
    fig = plt.figure(figsize=(11, 8.5))

    # ax is the axis on which the log2s will be plotted
    # if Transcript provided, igv is the axis for the IGV-like figure
    # if conifer data is present, cnf is the axis for that plot

    if 'conifer' in df_gene.columns and transcript is not None:
        # apportion subplots
        gs = GridSpec(3, 1, height_ratios=[0.475, 0.475, 0.05])
        ax = fig.add_subplot(gs[0])
        cnf = fig.add_subplot(gs[1], sharex=ax)
        igv = fig.add_subplot(gs[2], sharex=ax)
        igv.set_ylim((0, 1))
        igv.axis('off')
        subplots = {'log2' : ax, 'conifer' : cnf}
        # label axes
        cnf.set_xlabel('Position', fontsize=12)
        cnf.set_ylabel('CoNIFER Ratio', fontsize=12)
        ax.set_ylabel('Adjusted Mean of Log Ratio', fontsize=12)

    elif 'conifer' in df_gene.columns:
        # apportion subplots
        gs = GridSpec(2, 1, height_ratios=[0.5, 0.5])
        ax = fig.add_subplot(gs[0])
        cnf = fig.add_subplot(gs[1], sharex=ax)
        subplots = {'log2' : ax, 'conifer' : cnf}
        # label axes
        cnf.set_xlabel('Position', fontsize=12)
        cnf.set_ylabel('CoNIFER Ratio', fontsize=12)
        ax.set_ylabel('Adjusted Mean of Log Ratio', fontsize=12)

    elif transcript is not None:
        # apportion subplots
        gs = GridSpec(2, 1, height_ratios=[0.95, 0.05])
        ax = fig.add_subplot(gs[0])
        igv = fig.add_subplot(gs[1], sharex=ax)
        igv.set_ylim((0, 1))
        igv.axis('off')
        subplots = {'log2' : ax}
        # label axes
        ax.set_xlabel('Position', fontsize=12)
        ax.set_ylabel('Adjusted Mean of Log Ratio', fontsize=12)
    
    else:
        # apportion subplots
        ax = fig.add_subplot(1,1,1)
        subplots = {'log2' : ax}
        # label axes
        ax.set_xlabel('Position', fontsize=12)
        ax.set_ylabel('Adjusted Mean of Log Ratio', fontsize=12)

    # create CNV subplots
    exon_labels = {}
    for data_column, axis in subplots.items():
        df_subplot = df_gene[['mean_pos', data_column, 'exon']]
        df_subplot = df_subplot[df_subplot[data_column].notnull()]
        exon_labels[data_column] = create_gene_subplot(axis, df_subplot, data_column)

    # create the igv subplot
    if transcript is not None:
        create_transcript_subplot(igv, transcript)

    # set title
    title_parts = df_gene.iloc[0][['chr', 'gene', 'transcript']]
    ax.set_title("Chromosome {}: {}:{}".format(*title_parts))
    
    # save first figure
    fig.tight_layout()
    pdf.savefig()

    # if any log2s were cutoff by [-2,2] plot, add second plot covering full range
    rescaled = False

    for data_column, axis in subplots.items():
        min_log = df_gene[data_column].min()
        max_log = df_gene[data_column].max()

        # check if rescaling needed
        if min_log < -1 * y_scale or max_log > y_scale:
            axis.set_ylim((min_log - 0.1, max_log + 0.1))
            rescaled = True
            # adjust exon label y-coordinates
            y_mid = (max_log + min_log) / 2
            stretch_factor = (max_log - min_log) / (y_scale * 2)

            for label in exon_labels[data_column]:
                # extract info from old label and remove it
                (old_x, old_y) = label.get_position()
                exon = label.get_text()
                color = label.get_color()
                new_y = old_y * stretch_factor + y_mid
                label.remove()
                # create a new label at the new coordinates
                axis.text(old_x, new_y, exon, fontsize=10, color=color)            

    if rescaled:
        ax.set_title("Chromosome {}: {}:{} (Plot 2)".format(*title_parts))

        # save second figure
        pdf.savefig()

    plt.close()

def plot_main(pdf, df, title, min_log_ratio, rolling_window_size):
    """Saves to the pdf a plot of every point across all chromosomes present in df, flagging genes above min_log_ratio"""
    fig = plt.figure(figsize=(11, 8.5))

    # ax is the axis on which the log2s will be plotted
    # if conifer data is present, cnf is the axis for that plot
    
    if 'conifer' in df.columns:
        ax = fig.add_subplot(211)
        cnf = fig.add_subplot(212, sharex=ax)
        plots = {'log2' : ax, 'conifer' : cnf}
        # label axes
        cnf.set_xlabel('Chromosome', fontsize=12)
        cnf.set_ylabel('CoNIFER Ratio', fontsize=12)
        ax.set_ylabel('Adjusted Mean of Log Ratio', fontsize=12)

    else:
        ax = fig.add_subplot(1,1,1)
        plots = {'log2' : ax}
        # label axes
        ax.set_xlabel('Chromosome', fontsize=12)
        ax.set_ylabel('Adjusted Mean of Log Ratio', fontsize=12)

    # plot entries left to right from chrom 1 to chrom X
    chromosomes = natsorted(df['chr'].unique())

    x_tick_values = []  # used create axis labels for each chromosome
    v_line_coords = []  # used to draw lines between each chromosome

    # create a scatter plot of CNV ratios vs index for each chromosome
    for chrom in chromosomes:
        df_chrom = df[df['chr'] == chrom]
        indices = df_chrom.index.values
        x_tick_values.append(np.median(indices))
        v_line_coords.append(np.max(indices))

        for data_column, axis in plots.items():
            axis.scatter(indices, df_chrom[data_column], label=chrom, marker=',', s=0.2)

    # flag genes that are above or below log ratio threshold
    moved_labels = {}
    for data_column, axis in plots.items():
        flagged_df = df[[data_column, 'gene']][df[data_column].notnull()]
        flagged_genes = flag_genes(flagged_df, min_log_ratio, rolling_window_size)

        # add labels for flagged genes
        moved_labels[data_column] = []
        for gene, coord in flagged_genes.items():
            label_x = coord[0]
            # adjust text coordinates if it's off the scale
            if coord[1] < -0.98 * y_scale:
                label_y = -0.95 * y_scale
                moved = True
            elif coord[1] > 0.98 * y_scale:
                label_y = 0.95 * y_scale
                moved = True
            else:
                label_y = coord[1]
                moved = False

            label = axis.text(label_x, label_y, gene, fontsize=8, va='center', ha='left')
        
            if moved:
                moved_labels[data_column].append((label, coord[1]))
    
    for axis in plots.values():
        # set plot limits, ticks, and tick labels
        axis.set_ylim((-1 * y_scale, y_scale))
        axis.locator_params(axis='both', nbins=4)
        axis.tick_params(right=True, top=True, labelsize=10)
        axis.set_xticks(x_tick_values)
        axis.set_xticklabels(chromosomes, fontsize=8)

        # add vertical lines (skip line after last chrom)
        for coord in v_line_coords[0:-1]:
            axis.axvline(x=coord, alpha=0.1, color='black', linewidth=1)

        # add horizontal lines
        axis.grid(axis='y', which='major', linestyle='solid', linewidth=1)
        for coord in [min_log_ratio * -1, min_log_ratio]:
            axis.axhline(y=coord, alpha=0.5, color='black', linewidth=1, linestyle='dotted')

    # add title
    ax.set_title(title)

    # create legend
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles,
               labels,
               fontsize=10, 
               markerfirst=False, 
               markerscale=20, 
               borderaxespad=0.0,
               loc='center right')

    # use a tight layout to maximize area for plots
    fig.tight_layout()

    # shift the subplots so the legend doesn't overlap
    plt.subplots_adjust(right=0.92)

    # save main figure
    pdf.savefig()

    # if log2s were cutoff by [-2,2] plot, add second plot covering full range
    rescaled = False

    for data_column, axis in plots.items():
        min_log = df[data_column].min()
        max_log = df[data_column].max()
        if min_log < -1 * y_scale or max_log > y_scale:
            axis.set_ylim((min_log - 0.2, max_log + 0.2))
            rescaled = True

    if rescaled:
        for data_column, axis in plots.items():
            for label, new_y in moved_labels[data_column]:
                gene = label.get_text()
                label_coords = label.get_position()
                label.remove()
                axis.text(label_coords[0], new_y, gene, fontsize=8, va='center', ha='left')

        ax.set_title(title + ' (Plot 2)')
        pdf.savefig()

    plt.close()

def action(args):
    # load the data
    df = load_cnv_data(args.cnv_data)

    # if refgene is supplied, create a mapping of gene name to Transcript
    transcripts = {}
    if args.refgene:
        df_ref = UCSCTable(args.refgene).data
        rows = df_ref.to_dict(orient='records')
        transcripts.update({row['name2'] : Transcript(row) for row in rows})

    # infer title for main plot if needed
    if not args.title:
        args.title = extract_plot_title(args.cnv_data)

    # create plots and save to pdf-formatted outfile
    with PdfPages(args.outfile) as pdf:

        # create the main plot
        plot_main(pdf, df, args.title, args.min_log_ratio, args.window_size)

        # create a plot for each gene
        for gene in df['gene'].unique():
            # don't create a plot for 'intergenic'
            if gene == 'intergenic':
                continue
            gene_df = df[df['gene']==gene]
            # if refgene has an entry for this gene, pass along the Transcript to make the IGV plot
            if transcripts.has_key(gene):
                transcript = transcripts[gene]
            else:
                transcript = None

            plot_gene(pdf, gene_df, transcript=transcript)