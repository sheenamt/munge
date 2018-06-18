"""Calculate quality control metrics for UMI tags and consensus generation.
"""
import collections
import math
import os
import sys
import argparse 

import numpy as np
import pysam
import csv
import pandas as pd

#from bokeh.io import show, save 
from bokeh.plotting import figure,save, output_file
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
    gb = df.groupby(df.columns.tolist()).size().reset_index().rename(columns={0:'# reads containing UMI'})
#    print('total reads:', gb['# reads containing UMI'].sum())
#    print('umis:', gb['umi'].count())


    #Histogram 
    # prepare data for plotting
    counts = gb.groupby(['# reads containing UMI']).size().reset_index().rename(columns={0:'count'})
 #   print('double check the total:', counts['count'].sum())
#    x = [1, 2, 3, 4, 5, 6, 7, 8, 10, 12]
    read_nums = counts['# reads containing UMI'].tolist()
    umis = counts['count'].tolist()  # source=source, legend="count",

#    top =[138881, 38033, 4130, 2831, 123, 353, 5, 57, 2, 1]

#    source = ColumnDataSource(counts)

    #creates plot
    p = figure(plot_width=600, plot_height=600, title="Reads counts grouped by UMI counts ") 

    #turn off scientific notation, y ranges from 0 -> millions
    p.left[0].formatter.use_scientific = False

    # add renderer
    p.vbar(x=read_nums, bottom=0, width=0.5,top=umis, color="navy")
    #add hovering for read counts
    p.add_tools(HoverTool(tooltips=[("#Reads", "@top")]))

    #where to generate output

    output_file=('umi_metrics.html') #, mode='inline')

    # show or save results
    save(p)

