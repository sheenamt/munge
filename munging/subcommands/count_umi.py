"""Calculate quality control metrics for UMI tags and consensus generation.
"""
import pysam
import pandas as pd

import plotly
import plotly.graph_objs as go

def build_parser(parser):
    parser.add_argument('bam',
                        help='Path to bam')
    parser.add_argument('outfile',
                        help='Path to outfile')

def _get_umi_tag(rec):
    """Handle UMI and duplex tag retrieval.
    """
    for tag in ["RX", "XC"]:
        try:
            return rec.get_tag(tag)
        except KeyError:
            sys.exit("UMI tag not found")

def action(args):
    umi_bam = args.bam

    #create data frame of chrom:start:UMI from bam
    df = pd.DataFrame(({'umi':_get_umi_tag(rec), 'chrom': rec.reference_name, 'start': rec.reference_start, 'end': rec.reference_end} for rec in pysam.AlignmentFile(umi_bam,'rb', check_sq=False)))

    #group by chrom:start:stop:UMI and count duplicates, adding that as a column and reseting the index now that duplicates are removed
    gb = df.groupby(df.columns.tolist()).size().reset_index().rename(columns={0:'num_umi_pos'})

    #group by num_umi_pos to count how many reads were supported by various num_umi_pos
    counts = gb.groupby(['num_umi_pos']).size().reset_index().rename(columns={0:'count'})
    
    #setup bar chart
    barchart = go.Bar(x=counts['num_umi_pos'], 
                      y=counts['count'],
                      xaxis='x1',
                      yaxis='y1',
                      name='Reads'
                  )

    #setup table
    table = go.Table(
        header = dict(values=['# reads per UMI', 'Count'],
                      line = dict(color='#7D7F80'),
                      fill = dict(color='#a1c3d1'),
                      align = ['left'] * 5),
        cells = dict(values=[counts['num_umi_pos'], counts['count']],
                     line = dict(color='#7D7F80'),
                     fill = dict(color='#EDFAFF'),
                     align = ['left'] * 5),
        domain=dict(x=[.7, 1],
                    y=[0, 1]))             

    #adjust layout so barchart and table are next to each other
    layout = dict(xaxis1=dict( dict(domain=[0, .7], anchor='x1')),
                  yaxis1=dict( dict(domain=[0, 1], anchor='y1')),
                  title='Reads counts grouped by UMI counts')

    fig = go.Figure(data = [barchart,table], layout = layout)

    #plot
    plotly.offline.plot(fig, filename=args.outfile, auto_open=False)


