"""
find some insertion breakpoints that are lacked by Delly
"""

from __future__ import division
import sys
import pysam
import numpy as np
import re
import os
from sklearn.cluster import DBSCAN
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

min_clipped_len = 30

def check_clipped(read):
    flag = False
    start = True
    locus_estimate = 0
    clipped_len = 0
    for ci in read.cigar: # record left clipped length
        if ci[0] == 4 or ci[0] == 5 or ci[0] == 1:   # clipped or insertion
            clipped_len += int(ci[1])
        elif start and ci[0] == 0: # the locus start from xxM, should + length
            locus_estimate = int(ci[1])
        start = False
    if clipped_len >= min_clipped_len:
        flag = True
    return flag, locus_estimate


def ins_bps(bamname):
    position_dict = {}
    # position_list = []
    bamfile = pysam.AlignmentFile(filename = bamname, mode = 'rb')
    for read in bamfile:
        if read.mapping_quality < 20:
            continue
        # if read.reference_name != "NODE_23_length_642_cov_33.653097":
        #     continue
        # if read.is_unmapped == False and read.mate_is_unmapped == True: 
        if read.is_unmapped == False : 
            flag, locus_estimate = check_clipped(read)
            if flag:
                if read.reference_name not in position_dict:
                    position_dict[read.reference_name] = []
                position_dict[read.reference_name].append(read.reference_start + locus_estimate)
                # if read.reference_name == "NODE_23_length_642_cov_33.653097":
                #     print (read.reference_start, read.cigar, locus_estimate)
                # print (read)
    cluster_central = defaultdict(list)


    for chrom_name in position_dict: # for each chrom
        position_list = position_dict[chrom_name]
        # print (chrom_name, sorted(position_list))
        position_list = np.array(position_list).reshape(-1, 1)
        clustering = DBSCAN(eps=10, min_samples=2).fit(position_list)
        cluster_num = max(clustering.labels_) + 1

        if cluster_num == 0:
            continue
        
        save_clusters = []
        
        for i in range(cluster_num):
            save_clusters.append([])
        # print (cluster_num, len(position_list), len(clustering.labels_), position_list)
        for j in range(len(position_list)):
            save_clusters[clustering.labels_[j]].append(position_list[j])
        for i in range(cluster_num):
            central = round(np.median(save_clusters[i]))
            # print ("breakpoint cluster", chrom_name, central)
            # cluster_central.append([chrom_name, central])
            if central > 50  and contig_lengths[chrom_name] - central > 50:
                cluster_central[chrom_name].append(central)
        # print (position_list, clustering.labels_)
    # 
    for chrom_name in cluster_central:
        cluster_central[chrom_name] = sorted(cluster_central[chrom_name])
    print (cluster_central)
    return cluster_central

def get_contig_lengths(fasta_file):
    contig_lengths = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        contig_lengths[record.id] = len(record.seq)
    return contig_lengths

def split_contigs(cluster_central):
    output_handle = open(splited_contig, "w")
    record_dict = SeqIO.to_dict(SeqIO.parse(raw_contig, "fasta"))

    for seg_name in record_dict:
        record = record_dict[seg_name]
        if seg_name not in cluster_central:
            SeqIO.write(record, output_handle, "fasta")
        else:
            breaks = cluster_central[seg_name]
            intervals = []
            start = 0
            for point in breaks:
                intervals.append([start, point])
                start = point
            intervals.append([start, contig_lengths[seg_name]])

            for interval in intervals:
                split_seg = str(record.seq)[interval[0]:interval[1]]
                split_record = SeqRecord(Seq(split_seg), f"{seg_name}_{interval[0]}_{interval[1]}", '', '')
                SeqIO.write(split_record, output_handle, "fasta")

    output_handle.close()


if __name__ == "__main__":            
    bamname = sys.argv[1]
    raw_contig = sys.argv[2]
    splited_contig = sys.argv[3]

    contig_lengths = get_contig_lengths(raw_contig)
    cluster_central = ins_bps(bamname)
    split_contigs(cluster_central)