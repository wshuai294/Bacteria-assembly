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


def check_clipped(read):
    flag = False
    clipped_len = 0
    for ci in read.cigar: # record left clipped length
        if ci[0] == 4 or ci[0] == 5:
            clipped_len += int(ci[1])
    if clipped_len > 40:
        flag = True
    return flag


def calCrossReads(bamfile):
    dict_Interact_Big = {}
    for read in bamfile:
        if read.mapping_quality < 20:
            continue
        if read.is_unmapped == False and read.mate_is_unmapped == False \
            and read.flag < 2048: 
            if check_clipped(read):
                
                if read.reference_name == "NC_000913.3" \
                    and read.reference_start >  2947743 \
                    and read.reference_start <  2987249:
                    print (read.reference_name, read.reference_start)

            
bamname = sys.argv[1]
bamfile = pysam.AlignmentFile(filename = bamname, mode = 'rb')
calCrossReads(bamfile)