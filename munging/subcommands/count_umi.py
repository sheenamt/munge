"""Calculate quality control metrics for UMI tags and consensus generation.
"""
import pysam
import pandas as pd

from bokeh.plotting import figure,save, output_file, reset_output
from bokeh.models import ColumnDataSource,HoverTool

def build_parser(parser):
    parser.add_argument('bam',
                        help='Path to bam')

def _get_umi_tag(rec):
    """Handle UMI and duplex tag retrieval.
    """
    for tag in ["RX", "XC"]:
        try:
            return rec.get_tag(tag)
        except KeyError:
            pass

def action(args):
    umi_bam = args.bam
    
    #create data frame of chrom:start:UMI from bam
    df = pd.DataFrame(({'umi':_get_umi_tag(rec), 'chrom': rec.reference_name, 'start': rec.reference_start} for rec in pysam.AlignmentFile(umi_bam,'rb', check_sq=False)))

    #group by chrom:pos:start and count duplicates, adding that as a column and reseting the index now that duplicates are removed
    gb = df.groupby(df.columns.tolist()).size().reset_index().rename(columns={0:'num_umi_pos'})

    #group by num_umi_pos to count how many reads were supported by various num_umi_pos
    counts = gb.groupby(['num_umi_pos']).size().reset_index().rename(columns={0:'count'})

    #setup data for plotting
    source = ColumnDataSource(counts)

    #creates plot
    p = figure(plot_width=600, plot_height=600, title="Reads counts grouped by UMI counts") #, tools=[hover]) 

    #turn off scientific notation, y ranges from 0 -> millions
    p.left[0].formatter.use_scientific = False

    #bar plot
    p.vbar(x='num_umi_pos', top='count',bottom=0, width=0.5, color="green", source=source)
 
    #add hovering for read counts
    p.add_tools(HoverTool(tooltips=[("#Reads", "@count")]))

    #where to generate output and what to call the plot
    output_file('umi_metrics.html', title='UMI Metrics')

    #save results
    save(p)
