"""
Munges AnnotSV annotation of GRIDSS output
"""
import logging
import argparse
import pandas as pd
import sys
from munging.annotation import multi_split, chromosomes

log = logging.getLogger(__name__)
pd.options.display.width = 1000
pd.options.display.max_columns=30


def build_parser(parser):
    parser.add_argument('annotsv', 
                        help='A required input file')
    parser.add_argument('-o', '--outfile',
                        help='Output file', default=sys.stdout,
                        type=argparse.FileType('w'))

                                                                                
def parse_sv_event1(data):
    ''' Combine fields to make Event1 
    '''
    try:
        data['Event1']='chr'+str(chromosomes[data['SV chrom']])+':'+str(data['SV start'])
    except KeyError:
        print 'here:',data['SV chrom']

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
        try:
            chrom=a.split(':')
            data['Event2']='chr'+str(chromosomes[a[0]])+str(a[1:])
            data['Seq']=b
        except KeyError:
            print 'a', a
    else:
        try:
            chrom=b.split(':')
            data['Event2']='chr'+str(chromosomes[b[0]])+str(b[1:])
            data['Seq']=a
        except KeyError:
            print 'b', b
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
    data['Gene']=";".join(filter(None, [gene,promoter]))

    return pd.Series(data)

def parse_dgv(data):
    '''Combine DGV gain/loss columns into
    gain_n/gain_n_samples loss_n/loss_n_samples'''

    data['DGV_GAIN_found|tested']=str(data['DGV_GAIN_n_samples_with_SV'])+'|'+str(data['DGV_GAIN_n_samples_tested'])
    data['DGV_LOSS_found|tested']=str(data['DGV_LOSS_n_samples_with_SV'])+'|'+str(data['DGV_LOSS_n_samples_tested'])
    return pd.Series(data)

def parse_location(data):
    '''Split the annotsv location into two, if both are the same
    intron37-intron37 returns intron37
    intron37-intron39 returns intron37-intron39
    '''

    if data['location'] != '':
        a,b=data['location'].split('-')
        if a == b:
            data['location']=a
    return pd.Series(data)

def parse_repeats(data):
    ''' Combine left and right repeat info into one column'''

    repeat_left,repeat_right=None,None
    if data['Repeats_type_left'] != '':
        repeat_left=data['Repeats_type_left']+'[left]'
    if data['Repeats_type_right'] != '':
        repeat_right=data['Repeats_type_right']+'[right]'

    data['Repeats']=';'.join(filter(None,[repeat_left,repeat_right]))
    return pd.Series(data)

def parse_singleton(event):
    return [event['Event1'], event['Event2'], event['Gene'], event['Gene name'], event['location'],'SINGLETON EVENT', event['NM'], event['QUAL'], 'SINGLETON EVENT;'+event['FILTER'], event['1000g_event'],event['1000g_max_AF'],event['Repeats'],'SINGLETON EVENT',event['DGV_GAIN_found|tested'],event['DGV_LOSS_found|tested']]

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

    # set o and h events
    for sub_event in sub_events:
        if sub_event.endswith('o'):
            o_event = sub_event
        elif sub_event.endswith('h'):
            h_event = sub_event
    
    # collapse event sides into one result
    o_dict = collapse_event(event_df.loc[(event_df['ID']==o_event)])
    h_dict = collapse_event(event_df.loc[(event_df['ID']==h_event)])
    
    if o_dict and not h_dict:
        print 'only 1 event found for {}, probably due to quality: {}'.format(sub_events[0], [x for x in event_df['QUAL']])
        return parse_singleton(o_dict)
    elif h_dict and not o_dict:
        print 'only 1 event found for {}, probably due to quality: {}'.format(sub_events[0], [x for x in event_df['QUAL']])
        return parse_singleton(h_dict)
       

    # combine sub_events into one event
    o_event1 = event_df.loc[event_df['ID']==o_event,'Event1'].iloc[0]
    o_event2 = event_df.loc[event_df['ID']==o_event,'Event2'].iloc[0]
    h_event1 = event_df.loc[event_df['ID']==h_event,'Event1'].iloc[0]
    h_event2 = event_df.loc[event_df['ID']==h_event,'Event2'].iloc[0]

    #double check we're labeling event1 and 2 correctly:
    try:
        assert(o_event2 == h_event1)
        assert(o_event1 == h_event2)
    except:
        print "Calls did not match for events o {}/h {}, expected: o1 {} == h2 {}; o2 {} == h1 {}".format(o_event, h_event, o_event1, h_event2, o_event2, h_event1)
        return ['Error', parse_singleton(o_dict), parse_singleton(h_dict)]

    #Great, set things to o_event
    event1=o_event1
    event2=o_event2

    #Remove duplicate gene entries or set to 'Intergenic' if no gene is present
    gene1=';'.join([x for x in set([x for x in o_dict['Gene'].split(';')])]) or 'Intergenic'
    gene2=';'.join([x for x in set([x for x in h_dict['Gene'].split(';')])]) or 'Intergenic'
    
    location1 = o_dict['location']
    location2 = h_dict['location']
    repeats1 = o_dict['Repeats']
    repeats2 = h_dict['Repeats']

    # Logic here to resolve  nm, qual, filter, 1000g_event, 1000g_max_AF, dgv_gain, dgv_loss
    for k in o_dict.keys():
        if k == 'Gene' or k == 'location' or k == 'Repeats':
            pass
        else:
            val = None
            if o_dict[k] == h_dict[k]:
                val = o_dict[k]
            else:
                val = ';'.join(filter(None,[str(h_dict[k]),str(o_dict[k])]))
            if k == '1000g_event':
                thousandg_event = val
            elif k == '1000g_max_AF':
                thousandg_max_AF = val
            elif k == 'NM':
                nm = val
            elif k == 'QUAL':
                qual = val
            elif k == 'FILTER':
                vcf_filter = val
            elif k == 'DGV_GAIN_found|tested':
                dgv_gain = val
            elif k == 'DGV_LOSS_found|tested':
                dgv_loss = val
            
    # return here
    return [event1, event2, gene1, gene2, location1, location2, nm, qual, vcf_filter, thousandg_event, thousandg_max_AF, repeats1, repeats2, dgv_gain, dgv_loss]

def collapse_event(sub_event_df):
    ''' Collapses set of o or h columns into one event'''
    result_dict = {}
    for key in sub_event_df.columns: #headers:
        if len(sub_event_df[key])>0:
            #do not join empty strings or 'nan'
            result_dict[key]=';'.join(str(x) for x in sub_event_df[key].unique() if str(x) != 'nan' and str(x) !='' )
    return result_dict


def action(args):
    #Setup columns for output
    var_cols = ['Event1', 'Event2', 'Gene1','Gene2','location1','location2','NM','QUAL','FILTER','1000g_event', '1000g_max_AF', 'Repeats1','Repeats2','DGV_GAIN_found|tested','DGV_LOSS_found|tested']

    #Make dataframe of annotsv annotation
    try:
        annotsv_df=pd.read_csv(args.annotsv, delimiter='\t', index_col=False, usecols=['SV chrom','SV start','SV end', 'ID', 'ALT','Gene name','NM','QUAL',
                                                                                       'FILTER','INFO','location','promoters','1000g_event', '1000g_max_AF', 
                                                                                       'Repeats_type_left', 'Repeats_type_right',
                                                                                       'DGV_GAIN_n_samples_with_SV','DGV_GAIN_n_samples_tested',
                                                                                       'DGV_LOSS_n_samples_with_SV','DGV_LOSS_n_samples_tested'])
    except ValueError:
        args.outfile.write('\t'.join(var_cols) + '\n')
        sys.exit()

    annotsv_df.fillna('', inplace=True)
        
    #filter all calls less than 200 quality
    annotsv_df=annotsv_df[annotsv_df['QUAL']>=200]
    if annotsv_df.empty:
        annotsv_df.to_csv(args.outfile, index=False, columns=var_cols,sep='\t')
        sys.exit()

    #filter calls that are not chr1-23,X,Y
    
    #Parse the parts we care about
    annotsv_df=annotsv_df.apply(parse_sv_event1, axis=1).apply(parse_sv_alt, axis=1).apply(parse_gene_promoter,axis=1).apply(parse_dgv, axis=1).apply(parse_repeats,axis=1).apply(parse_info, axis=1).apply(parse_location, axis=1)


    # get list of unique Event_ID
    events_list = annotsv_df['EventID'].unique()
    event_results_list = []
    for event_id in events_list:
        # get subset dataframe with just that event_id
        current_event_df = annotsv_df.loc[(annotsv_df['EventID']==event_id)]
        if len(current_event_df) == 0:
            print("something went wrong...")
            sys.exit(1)
        # feed current event df copy to new function
        event_result = smoosh_event_into_one_line(current_event_df.copy())
        #Sometimes the second event has a lower quality score that was filtered out, print that to a log and move on
        if event_result[0]=='Error':
            event_results_list.append(event_result[1])
            event_results_list.append(event_result[2])
        else:
            event_results_list.append(event_result)

    var_cols = ['Event1', 'Event2', 'Gene1','Gene2','location1','location2','NM','QUAL','FILTER','1000g_event', '1000g_max_AF', 'Repeats1','Repeats2','DGV_GAIN_found|tested','DGV_LOSS_found|tested']
    output_df=pd.DataFrame(event_results_list,columns=var_cols)
    output_df.to_csv(args.outfile, index=False, columns=var_cols,sep='\t')


    
