import logging

from __init__ import __version__

log = logging.getLogger(__name__)

def get_location(chr, start, stop, **kwargs):
    """
    Format chr_loc for output (chr:start-end)
    """
    if start == stop:
        chr_loc = 'chr%s:%s' % (chr, start)
    else:
        chr_loc = 'chr%s:%s-%s' % (chr, start, stop)
    return chr_loc

def multi_split(source, splitlist):
    """
    Function to split a string given a string of multiple split points
    """
    output=[]
    atsplit=True
    if source==None:
        return None
    else:
        for char in source:
            if char in splitlist:
                atsplit=True
            else:
                if atsplit:
                    output.append(char)
                    atsplit=False
                else:
                    output[-1]=output[-1]+char
    return output

def split_chr_loc(d):

    """
    Function to parse chr_loc(chr:start-end) into chrm, start, end
    """
    output = multi_split(d, 'chr:-')
    #if SNP, there is not end position so set end=start
    try:
        end=output[2]
    except IndexError:
        end=output[1]
    chrm = output[0]
    start = output[1]
    return chrm, start, end

def split_string_in_two(data):
    """
    Return info from one column in two columns
    """
    if not data:
        return '-1', '-1'
    else:
        try:
            output=data.split(',')
            freq=output[0]
            count=output[1]
        except IndexError:
            output=data.split('|')
            freq=output[0]
            count=output[1]
    return freq, count

def build_variant_id(data):
    """
    Construct a variant_id from dict `d` containing keys ....
    """
    d={}
    if data[0]=='Position':
        return 'link'
    else:
        d['chromosome'], d['start'], d['end'] = split_chr_loc(data[0])
        d['ref_base'], d['var_base'] =  data[1], data[2]
        return '{chromosome}_{start}_{end}_{ref_base}_{var_base}'.format(**d)


