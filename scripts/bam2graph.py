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
    bamfile = pysam.AlignmentFile(filename = bam_name, mode = 'r')
    mean, sdev, rlen = getInsertSize(bamfile)
    insert_size = int(mean + 2*sdev)
    rlen = int(rlen)
    print (insert_size, rlen)
    for read in bamfile:
        if read.mapping_quality < 20:
            continue
        # if read.is_unmapped == False and read.mate_is_unmapped == False \
        #     and read.flag < 2048: 
            # if args['n'] == 1:
            #     read.reference_start = int(read.reference_name.split(':')[1].split('-')[0]) + read.reference_start
            #     read.next_reference_start = int(read.next_reference_name.split(':')[1].split('-')[0]) + read.next_reference_start
            # if (read.reference_name.split(':')[0] != read.next_reference_name.split(':')[0]):
            #     pass
            # else:
            #     if abs(read.next_reference_start - read.reference_start) > 3*insert_size:
            #         pass
            #     elif check_clipped(read):
            #         pass
            #     else:
            #         continue
        if (read.reference_name != read.next_reference_name):
            print (read.reference_name, read.reference_start, read.next_reference_name, read.next_reference_start)

bam_name = sys.argv[1]
calCrossReads(bam_name)