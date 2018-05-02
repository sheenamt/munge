"""Calculate quality control metrics for UMI tags and consensus generation.
"""
import collections
import math
import os
import sys
import argparse 

import numpy as np
import pysam
import yaml
import csv

def build_parser(parser):
    parser.add_argument('bam',
                        help='Path to bam')
    parser.add_argument('idxstats',
                        help='Path to idxstats of bam, from samtools')
    parser.add_argument('-o', '--outfile',
                        help='Output file', default=sys.stdout,
                        type=argparse.FileType('w'))

#Pulls from config
#get_sample_name 
#get_align_bam 

def _get_umi_tag(rec):
    """Handle UMI and duplex tag retrieval.
    """
    for tag in ["RX", "XC"]:
        try:
            return rec.get_tag(tag)
        except KeyError:
            pass

def safe_makedir(dname):
    """Make a directory if it doesn't exist, handling concurrent race conditions.
    """
    if not dname:
        return dname
    num_tries = 0
    max_tries = 5
    while not os.path.exists(dname):
        # we could get an error here if multiple processes are creating
        # the directory at the same time. Grr, concurrency.
        try:
            os.makedirs(dname)
        except OSError:
            if num_tries > max_tries:
                raise
            num_tries += 1
            time.sleep(2)
    return dname

def file_uptodate(fname, cmp_fname):
    """Check if a file exists, is non-empty and is more recent than cmp_fname.
    """
    try:
        return (file_exists(fname) and file_exists(cmp_fname) and
                os.path.getmtime(fname) >= os.path.getmtime(cmp_fname))
    except OSError:
        return False

def parse_idxstats(idxstats):
    """Return BAM index stats for the given file, using samtools idxstats.
    """
    AlignInfo = collections.namedtuple("AlignInfo", ["contig", "length", "aligned", "unaligned"])
    idxstats_out = open(idxstats, 'r')
    out = []
    for line in idxstats_out: #.split("\n"):
        if line.strip():
            contig, length, aligned, unaligned = line.split("\t")
            out.append(AlignInfo(contig, int(length), int(aligned), int(unaligned)))
    return out


def action(args):
#def run(_, data, out_dir):
    stats_file = args.outfile
    umi_bam = args.bam
    idxstats_file = args.idxstats
    out = {}
    total = 0
    mapped = 0
    duplicates = 0
    umi_reductions = []
    umi_counts = collections.defaultdict(int)
    with pysam.AlignmentFile(umi_bam, "rb", check_sq=False) as bam_iter:
        cur_counts = collections.defaultdict(int)
        cur_key = None
        for rec in bam_iter:
            total += 1
            umi = _get_umi_tag(rec)
            if umi and not rec.is_unmapped:
                mapped += 1
                if rec.is_duplicate:
                    duplicates += 1
                # chrom = bam_iter.getrname(rec.reference_id)
                # pos = rec.reference_start
                # key = (chrom, pos)
                key = umi
                #if this is a new chrm,pos
                if key != cur_key:
                    # update counts
                    if cur_counts:
                        for c in cur_counts.values():
                            umi_counts[c] += 1
                        total_seqs = sum(cur_counts.values())
                        umi_count = len(cur_counts)
                        umi_reductions.append(float(total_seqs) / umi_count)
                    # update current keys
                    cur_key = key
                    cur_counts = collections.defaultdict(int)
                cur_counts[umi] += 1
        if cur_counts:
            for c in cur_counts.values():
                umi_counts[c] += 1
            total_seqs = sum(cur_counts.values())
            umi_count = len(cur_counts)
            umi_reductions.append(float(total_seqs) / umi_count)
    out["umi_baseline_all"] = total
    out["umi_baseline_mapped"] = mapped
    out["umi_duplicates"] = duplicates
    out["umi_baseline_duplicate_pct"] = float(duplicates) / float(mapped) * 100.0
    out["umi_reduction_median"] = int(math.ceil(np.median(umi_reductions)))
    out["umi_reduction_max"] = int(max(umi_reductions))
    out["umi_counts"] = dict(umi_counts)
#    with open(stats_file, "w") as out_handle:

    # yaml.safe_dump({args.outfile: out}, out_handle,
    #                default_flow_style=False, allow_unicode=False)
    # return stats_file
    w = csv.writer(stats_file)
    for key, val in out.items():
        w.writerow([key, val])
