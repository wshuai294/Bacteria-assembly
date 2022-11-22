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

min_clipped_len = 30

def check_clipped(read):
    flag = False
    clipped_len = 0
    for ci in read.cigar: # record left clipped length
        if ci[0] == 4 or ci[0] == 5:
            clipped_len += int(ci[1])
    if clipped_len >= min_clipped_len:
        flag = True
    return flag


def ins_bps(bamname):
    position_dict = {}
    # position_list = []
    bamfile = pysam.AlignmentFile(filename = bamname, mode = 'rb')
    for read in bamfile:
        if read.mapping_quality < 20:
            continue
        if read.is_unmapped == False and read.mate_is_unmapped == True: 
            if check_clipped(read):
                if read.reference_name not in position_dict:
                    position_dict[read.reference_name] = []
                position_dict[read.reference_name].append(read.reference_start)
                # print (read)
    cluster_central = []


    for chrom_name in position_dict: # for each chrom
        position_list = position_dict[chrom_name]
        position_list = np.array(position_list).reshape(-1, 1)
        clustering = DBSCAN(eps=100, min_samples=2).fit(position_list)
        cluster_num = max(clustering.labels_) + 1
        
        save_clusters = []
        for i in range(cluster_num):
            save_clusters.append([])
        for j in range(len(position_list)):
            save_clusters[clustering.labels_[j]].append(position_list[j])
        for i in range(cluster_num):
            central = round(np.median(save_clusters[i]))
            print (chrom_name, central)
            cluster_central.append([chrom_name, central])
        # print (position_list, clustering.labels_)
    return cluster_central


if __name__ == "__main__":            
    bamname = sys.argv[1]
    ins_bps(bamname)