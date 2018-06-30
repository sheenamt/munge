"""Plot QC metrics, including on/off target coverage
"""
import pandas as pd
import sys
import os

from itertools import ifilter
from munging import filters,parsers
from munging.utils import walker, munge_pfx

import plotly
import plotly.graph_objs as go

def build_parser(parser):
    parser.add_argument('path',
                        help='Path to Quality_Analysis files')
    parser.add_argument('outfile',
                        help='Path to outfile')

def action(args):
    #grab the files
    filelist = ifilter(filters.quality_analysis, walker(args.path))

    #read in the data, adding the name into the df and skipping the 'version'
    df_list=[]
    for pfx_file in filelist:
        pfx = munge_pfx(pfx_file.fname)
        log_pfx=pfx['mini-pfx']
        df_list.append(pd.read_csv(os.path.join(pfx_file.dir, pfx_file.fname), sep='\t', skiprows=[2]).assign(SAMPLE=log_pfx))

    #concatenate them together
    big_df = pd.concat(df_list, ignore_index=True)
    #now, lets grab just the data we want 
    qc_df = big_df[['SAMPLE','MEAN_TARGET_COVERAGE','PCT_USABLE_BASES_ON_TARGET','PERCENT_DUPLICATION','PF_UNIQUE_READS']]

    data = [
        go.Bar(
            x=qc_df['SAMPLE'], # assign x as the dataframe column 'x'
            y=qc_df['PCT_USABLE_BASES_ON_TARGET']*qc_df['PF_UNIQUE_READS'],
            name='on target',
            xaxis='x1',
            yaxis='y1',

        ),
        go.Bar(
            x=qc_df['SAMPLE'],
            y=qc_df['PF_UNIQUE_READS']-qc_df['PCT_USABLE_BASES_ON_TARGET']*qc_df['PF_UNIQUE_READS'],
            name='off target',
            xaxis='x1',
            yaxis='y1',

        ),
        go.Bar(
            x=qc_df['SAMPLE'],
            y=qc_df['PF_UNIQUE_READS']/(1-qc_df['PERCENT_DUPLICATION'])-qc_df['PF_UNIQUE_READS'],
            name='duplicates',
            xaxis='x1',
            yaxis='y1',

        ),
    ]

    layout= {
        'title':'QC Metrics',
        'xaxis': {'type':'category'}, #required so mini sample names are strings instead of numbers
        'yaxis': {'hoverformat':',f'},
        'barmode':'stack'}
    
    fig = go.Figure(data=data, layout=layout)

    plotly.offline.plot(fig, filename=args.outfile, auto_open=False)
    # #create the scatter plot with blue circles
    # scatterplot = go.Scatter(
    #     x = msi_df.index,
    #     y = msi_df['msings_score'],
    #     name = 'msings_score',
    #     mode = 'markers',
    #     marker = dict(
    #         size = 10,
    #         color = '#a1c3d1',
    #         line = dict(
    #             width = 2,
    #             color = '#000000')
    #     ))

    # #scatter plot has a line at the threshold, which can be one or two values
    # min_line =  {'type': 'line',
    #               'x0': 0,
    #               'y0': min_thres,
    #               'x1': num_samples,
    #               'y1': min_thres,
    #               'line': {
    #                   'color': '#000000',
    #                   'width': 2,
    #                   'dash': 'dashdot'}
    #             }
    # if min_thres != max_thres:
    #     max_line =  {'type': 'line',
    #                  'x0': 0,
    #                  'y0': max_thres,
    #                  'x1': num_samples, 
    #                  'y1': max_thres,
    #                  'line': {
    #                      'color': '#000000',
    #                      'width': 2,
    #                      'dash': 'dashdot'}
    #             }
    #     lines = [min_line, max_line]
    # else:
    #     lines = [min_line]
    
    # #setup the layout so scatterplot and table are above/below  each other
    # layout= {
    #     'title':'MSI_STATUS',
    #     'xaxis': {'title': 'Samples',
    #               'range':[-.5,num_samples],
    #               'domain':[0, 1], #go across entire screen
    #               'anchor':'x1'},
    #     'yaxis': {'title': "msings score",
    #               'range': [0,1],
    #               'domain':[.7,1], #only take bottom portion of screen
    #               'anchor':'y1'},
    #     'shapes':lines}

    # #setup table
    # table = go.Table(

    #     header = dict(values = ['Sample_ID', 
    #                             'unstable_loci', 
    #                             'covered_loci', 
    #                             'msings_score', 
    #                             'msi_status', 
    #                             'tumor_mutation_burden'], #being explicit on purpose
    #                   line = dict(color='#7D7F80'),
    #                   fill = dict(color='#a1c3d1'),
    #                   align = ['left'] * 5),
    #     cells = dict(values=[msi_df.index,
    #                          msi_df.unstable_loci,
    #                          msi_df.covered_loci,
    #                          msi_df.msings_score,
    #                          msi_df.msi_status,
    #                          msi_df.tumor_mutation_burden],
    #                  line = dict(color='#7D7F80'),
    #                  fill = dict(color='#EDFAFF'),
    #                  align = ['left'] * 5),
    #     # domain=dict(x=[.6, 1], #side by side
    #     #             y=[0, 1]))             
    #     domain=dict(x=[0, 1], #above/belowe
    #                 y=[0, .5]))             

    # fig = go.Figure(data = [scatterplot,table], layout = layout)

    # #Create output file
    # plotly.offline.plot(fig, filename=args.outfile, auto_open=False)
