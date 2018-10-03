"""
Munges AnnotSV annotation of GRIDSS output
"""
import sys
import subprocess
import logging
import os
import argparse
import pandas as pd

from munging.annotation import multi_split

log = logging.getLogger(__name__)
pd.options.display.width = 180


def build_parser(parser):
    parser.add_argument('annotsv', 
                        help='A required input file')
    parser.add_argument('-o', '--outfile',
                        help='Output file', default=sys.stdout,
                        type=argparse.FileType('w'))


                                                                                
def parse_sv_event1(data):
    ''' Combine fields to make Event1 
    '''
    #Create Event1
    data['Event1']='chr'+str(data['SV chrom'])+':'+str(data['SV start'])
    return pd.Series(data)

def parse_sv_alt(data):
    ''' Split the ALT field from breakend format 
    ]7:55248960]A  to chr7:55248960
    A[7:55249011[  to chr7:55249011
    '''

    a,b=multi_split(data['ALT'],'[]')

    #Create Event2
    #the position has a : in it while the sequence does not
    if ':' in a:
        data['Event2']='chr'+a
        data['Seq']=b
    else:
        data['Event2']='chr'+b
        data['Seq']=a
    return pd.Series(data)

def parse_info(data):
    '''Get the EventID for each read '''
    #Parse read depth and SVtype
    info=dict(item.split('=') for item in data['INFO'].split(";") if "=" in item)
    data['EventID']=info['EVENT']
    return pd.Series(data)


def parse_gene_promoter(data):
    ''' Combine promoter and gene fields '''
    promoter, gene= None, None
    if data['promoters'] != '':
        promoter=str(data['promoters'])+'[Promoter]'
    if data['Gene name'] != '':
        gene=data['Gene name']

    if gene and promoter:
        data['Gene']=gene+';'+promoter
    elif gene:
        data['Gene']=gene
    elif promoter:
        data['Gene']=promoter
    return pd.Series(data)

def parse_dgv(data):
    '''Combine DGV gain/loss columns into
    gain_n/gain_n_samples loss_n/loss_n_samples'''
    data['DGV_GAIN_found|tested']=str(data['DGV_GAIN_n_samples_with_SV'])+'|'+str(data['DGV_GAIN_n_samples_tested'])
    data['DGV_LOSS_found|tested']=str(data['DGV_LOSS_n_samples_with_SV'])+'|'+str(data['DGV_LOSS_n_samples_tested'])
    return pd.Series(data)

def parse_repeats(data):
    ''' Combine left and right repeat info into one column'''
    repeat_left,repeat_right=None,None

    if data['Repeats_type_left'] != '':
        repeat_left=data['Repeats_type_left']+'[left]'
    if data['Repeats_type_right'] != '':
        repeat_right=data['Repeats_type_right']+'[right]'

    if repeat_left and repeat_right:
        data['Repeats']=';'.join([repeat_left,repeat_right])
    elif repeat_left:
        data['Repeats']=repeat_left
    elif repeat_right:
        data['Repeats']=repeat_right
    return pd.Series(data)

def smoosh_event_into_one_line(event_df):
    ''' Smooshes a multiline annotsv event into one line'''

    result = []
    event1,gene1,location1,repeats1 = None,None,None,None
    event2,gene2,location2,repeats2 = None,None,None,None
    nm = None
    qual = None
    vcf_filter = None
    thousandg_event = None
    thousandg_max_AF = None
    dgv_gain = None
    dgv_loss = None

    sub_events = event_df['ID'].unique()
    o_event = None
    h_event = None
    if len(sub_events) != 2:
        print("something went wrong later")
        sys.exit()
    # set o and h events
    for sub_event in sub_events:
        if sub_event.endswith('o'):
            o_event = sub_event
        elif sub_event.endswith('h'):
            h_event = sub_event

    # collapse event sides into one result
    print(event_df.loc[(event_df['ID']==o_event)])

    o_dict = collapse_event(event_df.loc[(event_df['ID']==o_event)])
    sys.exit()

    h_dict = collapse_event(event_df.loc[(event_df['ID']==h_event)])

    # combine sub_events into one event
    event1 = o_event
    event2 = h_event
    gene1 = o_dict['Gene']
    gene2 = h_dict['Gene']
    location1 = o_dict['location']
    location2 = h_dict['location']
    repeats1 = o_dict['Repeats']
    repeats2 = h_dict['Repeats']

    # Logic here to resolve  nm, qual, filter, 1000g_event, 1000g_max_AF, dgv_gain, dgv_loss
    for k in o_dict.keys():
        if k == 'gene' or k == 'location' or k == 'repeats':
            pass
        else:
            val = None
            if o_dict[k] == h_dict[k]:
                val = o_dict[k]
            else:
                val = str(o_dict[k]) + ";" + str(h_dict[k])
            if k == 'thousandg_event':
                thousandg_event = val
            if k == 'thousandg_max_AF':
                thousandg_max_AF = val
            if k == 'nm':
                nm = val
            if k == 'qual':
                qual = val
            if k == 'filter':
                filter = val
    
    # return here
    result = [event1, event2, gene1, gene2, location1, location2, nm, qual, filter, thousandg_event, thousandg_max_AF, repeats1, repeats2, dgv_gain, dgv_loss]


def collapse_event(sub_event_df):
    ''' Collapses set of o or h columns into one event'''
    result_dict = {}
    headers = ['Gene','NM','location','1000g_event', '1000g_max_AF', 'Repeats','DGV_GAIN_found/tested','DGV_LOSS_found/tested']
    result_dict['Gene'] = ";".join(sub_event_df['Gene'].unique())

    for key in headers:
        print(sub_event_df[key])
    #     result_dict[key]=';'.join([sub_event_df[key]]).unique()
    # print(result_dict['Gene'])
    return result_dict

def action(args):
    #Make dataframe of annotsv annotation
    annotsv_df=pd.read_csv(args.annotsv, delimiter='\t', index_col=False, usecols=['SV chrom','SV start','SV end', 'ID', 'ALT','Gene name','NM','QUAL',
                                                                                   'FILTER','INFO','location','promoters','1000g_event', '1000g_max_AF', 
                                                                                   'Repeats_type_left', 'Repeats_type_right',
                                                                                   'DGV_GAIN_n_samples_with_SV','DGV_GAIN_n_samples_tested',
                                                                                   'DGV_LOSS_n_samples_with_SV','DGV_LOSS_n_samples_tested'])
    annotsv_df.fillna('', inplace=True)

    #Parse the parts we care about
    #####
    #Should this part be done later, when we act on each 'EventID'?
    #####
    annotsv_df=annotsv_df.apply(parse_sv_event1, axis=1).apply(parse_sv_alt, axis=1).apply(parse_gene_promoter,axis=1).apply(parse_dgv, axis=1).apply(parse_repeats,axis=1).apply(parse_info, axis=1)

    #filter all calls less than 200 quality
    annotsv_df=annotsv_df[annotsv_df['QUAL']>=200]

    # get list of unique Event_ID
    events_list = annotsv_df['EventID'].unique()

    event_results_list = []
    for event_id in events_list:
        # get subset dataframe with just that event_id
        current_event_df = annotsv_df.loc[(annotsv_df['EventID']==event_id)]
        if len(current_event_df) == 0:
            print("something went wrong...")
            sys.exit()
        # feed current event df copy to new function
        event_result = smoosh_event_into_one_line(current_event_df.copy())
        event_results_list.append(event_result)

    # turn event_results_list into DF and/or output




    # for name,group in annotsv_output:
    #     print('name:',name)
    #     print('group:',group)

    print(annotsv_output)
    var_cals = ['Event1', 'Event2', 'ID','EventID','Gene','NM','QUAL','FILTER','location','1000g_event', '1000g_max_AF', 'Repeats','DGV_GAIN_found/tested','DGV_LOSS_found/tested']

 #   var_cals = ['Event1', 'Event2', 'ID','Gene','NM','QUAL','FILTER','location','1000g_event', '1000g_max_AF', 'Repeats','DGV_GAIN_found/tested','DGV_LOSS_found/tested']
#    annotsv_df.to_csv(args.outfile, index=False, columns=var_cals,sep='\t')
#    annotsv_output.to_csv(args.outfile, index=False, columns=var_cals,sep='\t')


    
