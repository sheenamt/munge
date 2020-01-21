"""Plot metrics for UMI tags and consensus generation.
"""
import pandas as pd

import plotly
import plotly.graph_objs as go

import sys

def build_parser(parser):
    parser.add_argument('groupumi_metrics',
                        help='Path to metric file from fgbio GroupReadsByUmi')
    parser.add_argument('out_html',
                        help='Path to html output')
    parser.add_argument('out_table',
                        help='Path to table output')

def action(args):
    umi_data = args.groupumi_metrics

    #create data frame of chrom:start:UMI from bam
    df = pd.read_csv(umi_data, sep='\t')
    
    #setup bar chart
    barchart = go.Bar(x=df['family_size'],
                      y=df['count'],
                      xaxis='x1',
                      yaxis='y1',
                      name='reads'
                  )

    #setup table
    table = go.Table(
        header = dict(values=['Family Size', 'Count'],
                      line = dict(color='#7D7F80'),
                      fill = dict(color='#a1c3d1'),
                      align = ['left'] * 5),
        cells = dict(values=[df['family_size'], df['count']],
                     line = dict(color='#7D7F80'),
                     fill = dict(color='#EDFAFF'),
                     align = ['left'] * 5),
        domain=dict(x=[.7, 1],
                    y=[0, 1]))             

    #adjust layout so barchart and table are next to each other
    layout= {'title':'Distribution of UMI family sizes observed during grouping',
             'xaxis': {'title': 'Family Size',
                       'domain':[0,.7], #dont' go across entire screen
                       'type':'category',
                       'anchor':'x1'},
             'yaxis': {'title': "Counts",
                       'domain':[0,1],
                       'hoverformat':',f',
                       'anchor':'y1'}}


    fig = go.Figure(data = [barchart,table], layout = layout)

    #plot
    plotly.offline.plot(fig, filename=args.out_html, auto_open=False)

    #    family_size Int The family size, or number of templates/read-pairs belonging to the family.
    #    count Count The number of families (or source molecules) observed with family_size observations.
    #    fraction Proportion The fraction of all families of all sizes that have this specific family_size.
    #    fraction_gt_or_eq_family_size Proportion The fraction of all families that have >= family_size

    #table
    df.head().to_csv(args.out_table, index=False, sep='\t')

