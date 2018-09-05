"""Plot QC metrics, including on/off target coverage
"""
#Standard library imports
import sys
import os
import re
from itertools import ifilter

#Libraries required to be installed 
from collections import namedtuple
import pandas as pd
import plotly
import plotly.graph_objs as go

from itertools import ifilter
from munging import filters,parsers
from munging.utils import walker, munge_pfx

def build_parser(parser):
    parser.add_argument('path',
                        help='Path to Quality_Analysis files')
    parser.add_argument('outfile',
                        help='Path to outfile')

def action(args):

    #read in the data, adding the name into the df and skipping the 'version'
    filelist = ifilter(filters.hs_file_finder, walker(args.path
    df_list=[]
    pd.set_option('display.width', 100)
    for pfx_file in filelist:
        pfx = munge_pfx(pfx_file.fname)
        log_pfx=pfx['mini-pfx']
        data=pd.read_csv(os.path.join(pfx_file.dir, pfx_file.fname), sep='\t',comment='#',error_bad_lines=False).assign(SAMPLE=log_pfx)
        
        df_list.append(data[0:1])

    #concatenate them together
    big_df = pd.concat(df_list, ignore_index=True)

    # #now, lets grab just the data we want 
    qc_df = big_df[['SAMPLE','MEAN_TARGET_COVERAGE',
                    #              'MEDIAN_TARGET_COVERAGE','MAX_TARGET_COVERAGE',
                    'PCT_USABLE_BASES_ON_TARGET',
                    'PF_UNIQUE_READS',]]


    #Setup the values we wish to plot
    qc_df['On Target Reads']=qc_df['PF_UNIQUE_READS']*qc_df['PCT_USABLE_BASES_ON_TARGET']
    qc_df['Off Target Reads']=qc_df['PF_UNIQUE_READS']-qc_df['On Target Reads']
    qc_df['On Target Reads']=qc_df['On Target Reads'].astype(int)
    qc_df['Off Target Reads']=qc_df['Off Target Reads'].astype(int)

    qc_df=qc_df.sort_values(by=['SAMPLE'])

    #Setup the plot
    data1 = go.Bar(
        x=qc_df['SAMPLE'], # assign x as the dataframe column 'x'
        y=qc_df['Off Target Reads'],
        name='Off Target Reads',
        xaxis='x1',
        yaxis='y1')
    data2=go.Bar(
        x=qc_df['SAMPLE'],
        y=qc_df['On Target Reads'],
        name='On Target Reads',
        xaxis='x1',
        yaxis='y1')
    layout= {
        'title':'QC Metrics',
        'xaxis': {'type':'category', #required so mini sample names are strings instead of numbers
                  'domain':[0, 1]},
        'yaxis': {'hoverformat':',f', #print real numbers, not 149.786k
                  'domain':[.7,1]}, #only take bottom portion of screen}, 
        'barmode':'stack'}
    
    #setup table
    table = go.Table(
        header = dict(values=['Sample ID', 'Mean Target Coverage','Total Read Pairs','On Target Reads', 'Off Target Reads'],
                      line = dict(color='#7D7F80'),
                      fill = dict(color='#a1c3d1'),
                      align = ['left'] * 5),
        cells = dict(values=[qc_df['SAMPLE'], qc_df['MEAN_TARGET_COVERAGE'],qc_df['PF_UNIQUE_READS'],qc_df['On Target Reads'],qc_df['Off Target Reads']],
                     line = dict(color='#7D7F80'),
                     fill = dict(color=['rgb(245,245,245)', #color for the first column, red if Target Coverage below 500
                                        ['rgba(0,250,0, 0.8)' if val>=500 else 'rgba(250,0,0, 0.8)' for val in qc_df['MEAN_TARGET_COVERAGE']] ]),
                     align = ['left'] * 5),
        domain=dict(x=[0, 1], #above/belowe
                    y=[0, .5]))             
    
    #Make the plot
    fig = go.Figure(data=[data1,data2,table], layout=layout)
    
    plotly.offline.plot(fig, filename=args.outfile, auto_open=false)
