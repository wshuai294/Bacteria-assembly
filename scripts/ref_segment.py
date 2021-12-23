import csv
import pysam
from pyfaidx import Fasta
import skbio
from skbio import DNA, TabularMSA
from skbio.alignment import local_pairwise_align_ssw
from skbio.alignment import StripedSmithWaterman
import numpy as np
import argparse
import sys
from sklearn.cluster import DBSCAN
import os


MIN_SEG_LEN = 50

class My_bkps():
    def __init__(self):
        self.all_bkp = []
        self.all_pos = {}
        self.cluster_bandwidth = 50
        self.ref_len = {}
        self.get_ref_len()

        self.segments = []
    
    def add(self, one_bkp):
        self.all_bkp.append(one_bkp)

    def add_pos(self, one_bkp):
        if one_bkp.from_ref in self.all_pos.keys():
            self.all_pos[one_bkp.from_ref].append(one_bkp.from_bkp)
        else:
            self.all_pos[one_bkp.from_ref]=[one_bkp.from_bkp]

        if one_bkp.to_ref in self.all_pos.keys():
            self.all_pos[one_bkp.to_ref].append(one_bkp.to_bkp)
        else:
            self.all_pos[one_bkp.to_ref]=[one_bkp.to_bkp]

    def read_bkp(self):       
        with open(acc_bkp_file, newline='') as csvfile:
            spamreader = csv.reader(csvfile, delimiter=',')
            for row in spamreader:
                if row[0] == "from_ref": #header line
                    continue
                one_bkp = Bkp_Record(row)
                print(', '.join(row))
                self.add(one_bkp)
                self.add_pos(one_bkp)

    def cluster_pos(self):
        for ref in self.all_pos:
            ref_pos_list = self.all_pos[ref]
            ref_rep_pos_list = self.dbscan_clu(ref_pos_list)
            self.all_pos[ref] = ref_rep_pos_list

    def dbscan_clu(self, ref_pos_list):
        ref_rep_pos_list = []
        XY = np.array(ref_pos_list).reshape(-1, 1)
        db = DBSCAN(eps = self.cluster_bandwidth, min_samples = 1).fit(XY)
        labels = db.labels_
        cluster_label_dict = {}
        lab = labels.tolist()
        for i in range(0,len(lab)):
            if lab[i] not in cluster_label_dict:
                ref_rep_pos_list.append(ref_pos_list[i])
                cluster_label_dict[lab[i]] = 1

        # print (ref_pos_list)
        # print (ref_rep_pos_list)
        return sorted(ref_rep_pos_list)

    def get_ref_len(self):
        f = open(ref_file+".fai", 'r')
        for line in f:
            array = line.split()   
            self.ref_len[array[0]] = int(array[1])          

    def get_segments(self):
        for ref in self.ref_len.keys():
            if ref not in self.all_pos:
                self.segments.append([ref, 1, self.ref_len[ref]])
            else:
                start = 1
                for pos in self.all_pos[ref]:
                    if pos - start < MIN_SEG_LEN:
                        continue
                    self.segments.append([ref, start, pos])
                    start = pos
                if self.ref_len[ref] - start >= MIN_SEG_LEN:
                    self.segments.append([ref, start, self.ref_len[ref]])
        # print (self.segments)
        # self.remove_unmapped_segs()
        self.get_segments_fasta()

    def get_bed_file(self):
        fout = open(bed_file, 'w')
        for seg in self.segments:
            print (f"{seg[0]}:{seg[1]}-{seg[2]}", file = fout)
        fout.close()

    def get_segments_fasta(self):
        self.get_bed_file()
        order = f"samtools faidx -r {bed_file} {ref_file} > {ref_seg_file}"
        os.system(order)

    def remove_unmapped_segs(self):
        for i in range(len(self.segments)):
            self.segments[i].append(0)

        for line in open(depth_file):
            array = line.strip().split()
            seg_name = array[0]
            pos = int(array[1])
            dp = int(array[2])
            for i in range(len(self.segments)):
                if self.segments[i][0] ==  seg_name and pos >= self.segments[i][1] and pos <= self.segments[i][2]:
                    self.segments[i][3] += 1
                    break
        for i in range(len(self.segments)):
            mapped_ratio = float(self.segments[i][3])/abs(self.segments[i][2]-self.segments[i][1])
            print (self.segments[i], mapped_ratio)



class Bkp_Record(object):
    def __init__(self, bkp_row):
        self.from_ref = bkp_row[0]
        self.from_bkp = int(bkp_row[1])
        self.to_ref = bkp_row[2]
        self.to_bkp = int(bkp_row[3])



            
        # writer = csv.writer(f)
        # header = ['from_ref','to_ref','from_pos','to_pos','from_side','to_side',\
        # 'if_reverse','read_seq','ref_seq','similarity']
        # writer.writerow(header)

# ref_file = "/mnt/d/breakpoints/assembly/simulation/ASM59784v1_genome.fa"
# ref_seg_file = "/mnt/d/breakpoints/assembly/simulation/ASM59784v1_genome.seg.fa"
# bed_file = "/mnt/d/breakpoints/assembly/simulation/ASM59784v1.bed"
# acc_bkp_file = "/mnt/d/breakpoints/assembly/simulation/ASM59784v1.acc.txt"


ref_file = sys.argv[1]
ref_seg_file = sys.argv[2]
acc_bkp_file = sys.argv[3]
depth_file = sys.argv[4]

bed_file = ref_file + ".bed"





my_bkps = My_bkps()
my_bkps.read_bkp()
my_bkps.cluster_pos()
my_bkps.get_segments()