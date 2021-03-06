"""
Munges AnnotSV annotation of GRIDSS output
"""
import logging
import argparse
import pandas as pd
import sys
import csv
from munging.annotation import multi_split, chromosomes
from intervaltree import Interval, IntervalTree

log = logging.getLogger(__name__)
pd.options.display.width = 1000
pd.options.display.max_columns=30


def build_parser(parser):
    parser.add_argument('annotsv', 
                        help='A required input file')
    parser.add_argument('-q','--quality_filter',type=int,
                        default=200,
                        help='Threshold for quality filter, 200 default')
    parser.add_argument('-s','--report_singletons',action='store_true',
                        help='Report singleton breakpoints. Results in many more output lines')
    parser.add_argument('-c','--capture_file',
                        help='Assay BED file. Allows annotations for on-target or off')
    parser.add_argument('-f','--flagged_fusions',
                        help='Assay flagged fusions file. Allows flagging of clinically-relevant gene fusions')                    
    parser.add_argument('-o', '--outfile',
                        help='Output file', default=sys.stdout,
                        type=argparse.FileType('w'))

                                                                                
def parse_length(event1,event2):
    '''Create leght if on same chrom,
    otherwise print CTX (inter-chromosomal translocation)
    '''
    chr1,e1=event1.split(':')
    chr2,e2=event2.split(':')

    #If on same chrome, print length
    if chr1==chr2:
        length=abs(int(e1)-int(e2))
    #Otherwise print CTX
    else:
        length='CTX'
    return length

def parse_quality(df, quality):
    ''' Remove calls that do not meet quality threshold
    '''
    ids=df.loc[df['QUAL']<=quality,'ID']
    #If all calls are above threshold, return the original df
    if ids.empty:
        return df
    else:
        ids_to_remove=[]
        for gid in ids:
            ids_to_remove.append(gid[:-1]+'h')
            ids_to_remove.append(gid[:-1]+'o')
        ids_to_remove=set(ids_to_remove)
        
        df=df[~df.ID.str.contains('|'.join(ids_to_remove))]
        return df

def parse_sv_event1(data):
    ''' Combine fields to make Event1 
    '''
    #Only process chr1-23, X, Y
    try:
        data['Event1']='chr'+str(chromosomes[data['SV chrom']])+':'+str(data['SV start'])
    except KeyError:
        pass

    return pd.Series(data)

def parse_sv_alt(data):
    ''' Split the ALT field from breakend format 
    ]7:55248960]A  to chr7:55248960
    A[7:55249011[  to chr7:55249011
    '''
    if data['ID'].endswith(('o','h')):
        a,b=multi_split(data['ALT'],'[]')
        #Create Event2
        #the position has a : in it while the sequence does not
        #Only process chr1-23, X, Y
        if ':' in a:
            try:
                chrom=a.split(':')
                data['Event2']='chr'+str(chromosomes[a[0]])+str(a[1:])
                data['Seq']=b
            except KeyError:
                pass
        else:
            try:
                chrom=b.split(':')
                data['Event2']='chr'+str(chromosomes[b[0]])+str(b[1:])
                data['Seq']=a
            except KeyError:
                pass
    elif data['ID'].endswith('b'):
        data['Event2']='SingleBreakEnd'
        data['Seq']=data['ALT']
    else:
        print("some weird data was encountered", data)
        sys.exit()
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

def split_genes(gene_field):
    """Undoes multi-gene concatenation and removes the 'Promoter' tag"""
    return str(gene_field).replace('[Promoter]','').split(';')

def parse_event_type(event):
    """Returns 'GENE_FUSION', 'INTRAGENE', 'INTERGENIC', or 'SINGLETON' depending on contents of event"""
    
    if event['Event2'] == 'SingleBreakEnd':
        return 'SINGLETON'

    # undo multi gene concatenation and remove 'Promoter' label
    genes_1 = split_genes(event['Gene1'])
    genes_2 = split_genes(event['Gene2'])

    # if there is no overlap between the set of genes for Gene1 and Gene2, label GENE_FUSION
    if not (set(genes_1) & set(genes_2)):
        return 'GENE_FUSION'
    # if both gene1 and gene2 are intergenic, label INTERGENIC
    elif genes_1 == ['Intergenic'] and genes_2 == ['Intergenic']:
        return 'INTERGENIC'
    # otherwise, label INTRAGENIC
    return 'INTRAGENIC'


def build_capture_trees(bed_file):
    """Reads BED file and returns a dict of 25 interval trees (1/per chromosome)"""
    # read the BED file
    bed_df = pd.read_csv(bed_file, usecols=[0,1,2,3], sep='\t',
        names=['chrom', 'start', 'stop', 'gene'],
        dtype={'chrom':str, 'start' : int, 'stop' : int, 'gene' : str})

    # Intervals exclude the end point, so increment all stops
    bed_df['stop'] += 1

    # build an interval tree for each chromosome
    trees = {}
    chroms = [str(i) for i in range(1,23)] + ['X', 'Y']
    for c in chroms:
       # select rows for c
       chrom_df = bed_df[bed_df.chrom == c][['start', 'stop', 'gene']]
       # convert rows to series of tuples
       chrom_tuples = [tuple(x) for x in chrom_df.values]
       # build the interval tree from tuples
       trees[c] = IntervalTree.from_tuples(chrom_tuples)
    
    return trees


def parse_capture_intent(event, interval_trees):
    """returns YES or NO depending on whether positions in event are found in the interval_trees"""
    
    # separate Event1 and Event2 into [<chr>, <pos>]
    pos1 = event['Event1'].replace('chr', '').split(':')
    pos2 = event['Event2'].replace('chr', '').split(':')

    # if pos1 is found in interval_trees label YES
    try:
         tree1 = interval_trees[pos1[0]]
         if len(tree1[int(pos1[1])]) > 0:
             return 'YES'
    except KeyError:
        pass

    # if pos2 is found in interval_trees label YES
    try:
         tree2 = interval_trees[pos2[0]]
         if len(tree2[int(pos2[1])]) > 0:
             return 'YES'
    except KeyError:
        pass

    # otherwise, label NO
    return 'NO'


def parse_flagged_fusions_file(flagged_fusions_file):
    """
    Uses the contents of the assay flagged fusions file to return a searchable data structure
    of annotated region pairs.

    Returns a dictionary of 24 IntervalTrees, one for each chromosome. Each IntervalTree
    is composed of Intervals whose boundaries are defined by [region1] and whose contents
    are themselves IntervalTrees.
    
    These secondary trees are composed of intervals whose boundaries are defined by [region2] and
    whose contents are [notes]
    """
    fusions_df = pd.read_csv(flagged_fusions_file, sep='\t', header=0, dtype={'chrom1':str, 'chrom2':str})

    # build an interval tree for each chromosome
    trees = {}
    chroms = ['chr' + str(i) for i in range(1,23)] + ['chrX', 'chrY']
    for c in chroms:
        
        # fiter the fusions_df on chrom1
        chrom_df = fusions_df[fusions_df['chrom1'] == c]

        trees[c] = IntervalTree()
        
        # create a collection of subtrees for each unique region1 on c
        regions = {}
        for _, row in chrom_df.iterrows():
            region1 = row['region1']

            # if this region hasn't been seen yet, create a dict for each chromosomal tree
            if region1 not in regions.keys():
                regions[region1] =  {}

            # the intervals of the subtrees are defined by region2
            chrom2 = row['chrom2']
            interval2 = row['region2'].split('-')
            start2 = int(interval2[0])
            stop2 = int(interval2[1]) + 1

            # if this is the first fusions partner on chrom2, create a new IntervalTree
            if chrom2 not in regions[region1].keys():
                regions[region1][chrom2] = IntervalTree()

            # add the interval for region2 
            regions[region1][chrom2][start2:stop2] = row['notes']
        
        # populate the primary interval tree for this chromosome
        for region1 in regions.keys():
            interval1 = region1.split('-')
            start1 = int(interval1[0])
            stop1 = int(interval1[1]) + 1
            trees[c][start1:stop1] = regions[region1]

    return trees


def parse_clinical_fusions(event, fusion_partners):
    """
    Returns a custom note if event is represented in the fusion_partners data structure
    Otherwise returns 'NO'.
    """
    # only consider gene fusion events
    if event['Type'] == 'GENE_FUSION':
        # find chrom/pos for each event; if they do not exist, it is not a flagged fusion
        try:
            region1 = event['Event1'].split(':')
            chrom1 = region1[0]
            pos1 = int(region1[1])
            region2 = event['Event2'].split(':')
            chrom2 = region2[0]
            pos2 = int(region2[1])
        except IndexError:
            return 'NO'

        # check if Event1 corresponds to a gene of interest and Event2 to one of its fusion partners
        # if so, return the pair's custom note
        for interval1 in fusion_partners[chrom1][pos1]:
            try:
                subtrees = interval1[2]
                for interval2 in subtrees[chrom2][pos2]:
                    return parse_flagged_fusion_note(interval2[2])

            except KeyError:
                continue

        # check if Event2 corresponds to a gene of interest and Event1 to one of its fusion partners
        # if so, return the pair's custom note
        for interval2 in fusion_partners[chrom2][pos2]:
            try:
                subtrees = interval2[2]
                for interval1 in subtrees[chrom1][pos1]:
                    return parse_flagged_fusion_note(interval1[2])
                    
            except KeyError:
                continue

    return 'NO'


def parse_flagged_fusion_note(note):
    """
    Modifies a flagged_fusion_note for display in a workbook
    """
    if 'Quiver' in note:
        fusion = note.split(':')[1].split('-')
        gene1 = fusion[0]
        gene2 = fusion[1]
        output = '=HYPERLINK("http://quiver.archerdx.com/results?query={0}%3A{1}", "{2}")'
        return output.format(gene1, gene2, note)
    
    return note

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
    b_event = None
    
    # set o,h and b events
    for sub_event in sub_events:
        if sub_event.endswith('o'):
            o_event = sub_event
        elif sub_event.endswith('h'):
            h_event = sub_event
        elif sub_event.endswith('b'):
            b_event = sub_event
        else:
            print("something went wrong...")
            print(sub_event)
            sys.exit()
    # collapse event sides into one result
    o_dict = collapse_event(event_df.loc[(event_df['ID']==o_event)])
    h_dict = collapse_event(event_df.loc[(event_df['ID']==h_event)])
    b_dict = collapse_event(event_df.loc[(event_df['ID']==b_event)])
    
    if o_dict and not h_dict:
        print 'only 1 event found for {}, probably due to location {}'.format(sub_events[0], list(set([x for x in event_df['ALT']])))
        return parse_singleton(o_dict)
    elif h_dict and not o_dict:
        print 'only 1 event found for {}, probably due to location {}'.format(sub_events[0], list(set([x for x in event_df['ALT']])))
        return parse_singleton(h_dict)
    elif b_dict:
        return parse_singleton(b_dict)
        
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

    #Calucate lenght, or print CTX if interchromosomal 
    length=parse_length(event1,event2)

    
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
    return [event1, event2, gene1, gene2, location1, location2, length, nm, qual, vcf_filter, thousandg_event, thousandg_max_AF, repeats1, repeats2, dgv_gain, dgv_loss]

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
                                                                                       'DGV_LOSS_n_samples_with_SV','DGV_LOSS_n_samples_tested',])
                                                                                       
    except ValueError:
        args.outfile.write('\t'.join(var_cols) + '\n')
        sys.exit()
    annotsv_df.fillna('', inplace=True)
    #filter all calls less than 200 quality
    annotsv_df=parse_quality(annotsv_df, quality=args.quality_filter)
    if annotsv_df.empty:
        annotsv_df.to_csv(args.outfile, index=False, columns=var_cols,sep='\t')
        sys.exit()
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
            #If the events were not a real chrom, it will be empty. Skip those
            if event_result[1][0] and event_result[2][0]:
                event_results_list.append(event_result[1])
                event_results_list.append(event_result[2])
        else:
            event_results_list.append(event_result)
    var_cols = ['Event1', 'Event2', 'Gene1','Gene2','location1','location2','Length', 'NM','QUAL','FILTER','1000g_event', '1000g_max_AF', 'Repeats1','Repeats2','DGV_GAIN_found|tested','DGV_LOSS_found|tested']
    output_df=pd.DataFrame(event_results_list,columns=var_cols)

    # add column for event 'Type'
    output_df['Type'] = output_df.apply(parse_event_type, axis=1)

    # if requested, add column for capture intent
    if args.capture_file:
        capture_trees = build_capture_trees(args.capture_file)
        output_df['Intended_For_Capture'] = output_df.apply(parse_capture_intent, interval_trees=capture_trees, axis=1)

    # if requested, add column for clinically relevant fusions
    if args.flagged_fusions:
        fp_trees = parse_flagged_fusions_file(args.flagged_fusions)
        output_df['Flagged_Fusions'] = output_df.apply(parse_clinical_fusions, fusion_partners=fp_trees, axis=1)

    # filter out singletons unless otherwise requested
    if not args.report_singletons:
        output_df=output_df[(output_df['Event2'] !='SingleBreakEnd') & (output_df['location2'] !='SINGLETON EVENT')]

    # save output to file
    output_df.to_csv(args.outfile, index=False, sep='\t', quoting=csv.QUOTE_NONE)

