from __future__ import division
import sys
import math
import pysam
import numpy as np
from sklearn.cluster import DBSCAN
import random
import re
import os
import multiprocessing
import argparse
import time
from get_raw_bkp import getInsertSize


def calCrossReads(bam_name):
    f = open(graph, "w")
    edge_dict = {}
    bamfile = pysam.AlignmentFile(filename = bam_name, mode = 'rb')
    mean, sdev, rlen = getInsertSize(bamfile)
    rlen = int(rlen)
    insert_size = int(mean + 2*sdev) + 2 * rlen
    
    for read in bamfile.fetch():
        if read.is_unmapped  or read.mate_is_unmapped:
            continue 
        if read.mapping_quality < min_q:
            continue
        if read.has_tag('XA'):
            continue
        if (read.reference_name == read.next_reference_name):
            continue

        ref_len = int(read.reference_name.split(':')[1].split('-')[1]) - int(read.reference_name.split(':')[1].split('-')[0])
        mate_len = int(read.next_reference_name.split(':')[1].split('-')[1]) - int(read.next_reference_name.split(':')[1].split('-')[0])

        if not (abs(read.reference_start) < insert_size or abs(ref_len - read.reference_start)< insert_size):
            continue
        if not (abs(read.next_reference_start) < insert_size or abs(mate_len - read.next_reference_start)< insert_size):
            continue

        if read.is_reverse:
            refname = read.reference_name + " -"
        else:
            refname = read.reference_name + " +"
        if not read.mate_is_reverse:
            mate_ref = read.next_reference_name + " -"
        else:
            mate_ref = read.next_reference_name + " +"
           

        edge = refname + " " + mate_ref

        if edge not in edge_dict:
            edge_dict[edge] = 1
        else:
            edge_dict[edge] += 1

    for edge in edge_dict:
        print (edge, edge_dict[edge], file = f)
        print (edge, edge_dict[edge])
    f.close()


min_q = 20
bam_name = sys.argv[1]
graph = sys.argv[2]
print ("start extract edges...")
calCrossReads(bam_name)