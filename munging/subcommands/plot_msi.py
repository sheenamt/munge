"""Plot all MSI outputs, including small table of scores
"""
import pandas as pd
import sys

import plotly
import plotly.graph_objs as go

def build_parser(parser):
    parser.add_argument('MSI_Analysis',
                        help='Path to combined MSI_Analysis file')
    parser.add_argument('outfile',
                        help='Path to outfile')
    parser.add_argument('-t','--msi_threshold',
                        nargs='+',
                        help='MSI score threshold or range of two thresholds, default of 0.2')
    

def parse_thresholds(threshold_list):
    """ Parse the user defined threshold list"""

    if len(threshold_list) == 1:
        min_thres=float(threshold_list[0])
        max_thres=float(threshold_list[0])
    elif len(threshold_list) == 2:
        min_thres=float(min(threshold_list))
        max_thres=float(max(threshold_list))
    else:
        raise ValueError("Wrong number of -t thresholds given")
    return min_thres, max_thres

def action(args):
    #read in the data
    msi_data = args.MSI_Analysis
    
    #setup the threshold for pos/neg
    if args.msi_threshold:
        threshold=args.msi_threshold
    else:
        threshold=[0.2,0.2]

    min_thres, max_thres = parse_thresholds(threshold)

    #create data frame of MSI data
    df = pd.read_csv(msi_data, sep='\t', index_col='Position')

    #transpose data so it is grouped by sample
    df_t=df.transpose()

    #Replace the tumor burden version, so that the script can be used for all versions
    df_t.columns = df_t.columns.str.replace('tumor_mutation_burden.*','tumor_mutation_burden')

    #grab all info besides site specific calls
    msi_df = df_t[df_t.columns[0:5]]

    #figure out number of samples, used in setting up plots
    num_samples = len(msi_df)

    #create the scatter plot with blue circles
    scatterplot = go.Scatter(
        x = msi_df.index,
        y = msi_df['msings_score'],
        name = 'msings_score',
        mode = 'markers',
        marker = dict(
            size = 10,
            color = '#a1c3d1',
            line = dict(
                width = 2,
                color = '#000000')
        ))

    #scatter plot has a line at the threshold, which can be one or two values
    min_line =  {'type': 'line',
                  'x0': 0,
                  'y0': min_thres,
                  'x1': num_samples,
                  'y1': min_thres,
                  'line': {
                      'color': '#000000',
                      'width': 2,
                      'dash': 'dashdot'}
                }
    if min_thres != max_thres:
        max_line =  {'type': 'line',
                     'x0': 0,
                     'y0': max_thres,
                     'x1': num_samples, 
                     'y1': max_thres,
                     'line': {
                         'color': '#000000',
                         'width': 2,
                         'dash': 'dashdot'}
                }
        lines = [min_line, max_line]
    else:
        lines = [min_line]
    
    #setup the layout so scatterplot and table are above/below  each other
    layout= {
        'title':'MSI_STATUS',
        'xaxis': {'title': 'Samples',
                  'range':[-.5,num_samples],
                  'domain':[0, 1], #go across entire screen
                  'type':'category',
                  'anchor':'x1'},
        'yaxis': {'title': "msings score",
                  'range': [0,1],
                  'domain':[.7,1], #only take bottom portion of screen
                  'anchor':'y1'},
        'shapes':lines}

    #setup table
    table = go.Table(
        header = dict(values = ['Sample_ID', 
                                'unstable_loci', 
                                'covered_loci', 
                                'msings_score', 
                                'msi_status', 
                                'tumor_mutation_burden'], #being explicit on purpose
                      line = dict(color='#7D7F80'),
                      fill = dict(color='#a1c3d1'),
                      align = ['left'] * 5),
        cells = dict(values=[msi_df.index,
                             msi_df.unstable_loci,
                             msi_df.covered_loci,
                             msi_df.msings_score,
                             msi_df.msi_status,
                             msi_df.tumor_mutation_burden],
                     line = dict(color='#7D7F80'),
                     fill = dict(color='#EDFAFF'),
                     align = ['left'] * 5),
        # domain=dict(x=[.6, 1], #side by side
        #             y=[0, 1]))             
        domain=dict(x=[0, 1], #above/belowe
                    y=[0, .5]))             

    fig = go.Figure(data = [scatterplot,table], layout = layout)

    #Create output file
    plotly.offline.plot(fig, filename=args.outfile, auto_open=False)
