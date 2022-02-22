import csv
import pysam
from pyfaidx import Fasta
import numpy as np
import argparse
import sys
from sklearn.cluster import DBSCAN
import os
from pysam import VariantFile
import re




def check_clipped(read):
    flag = False
    start = True
    clipped_len = 0
    pos_clip = 0
    for ci in read.cigar: # record left clipped length
        if ci[0] == 4 or ci[0] == 5:
            clipped_len += int(ci[1])
            start = False
        elif start == True:
            pos_clip += int(ci[1])
    if clipped_len > 50:
        flag = True
    return flag, pos_clip

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
                # print(', '.join(row))
                self.add(one_bkp)
                self.add_pos(one_bkp)

    def read_lumpy_vcf(self):
        # vcffile = "/mnt/d/breakpoints/assembly/simulation/assembly_test/test.sv.vcf"
        vcffile = acc_bkp_file
        if not os.path.isfile(vcffile):
            print ("no sv vcf.")
            return 0
        for line in open(vcffile, "r"):
            line = line.strip()
            if line[0] == "#":
                continue

            # """
            if re.search("IMPRECISE", line):
                continue


            len_r = re.search(";SVLEN=(.*?);", line)
            array = line.split()
            if len_r:
                sv_len = abs(int(len_r.group(1)))
                if sv_len < min_sv_len:
                    continue
                # print (sv_len)
                # """
                
                chrom = array[0]
                pos = int(array[1])
                self.add_lumpy(chrom, pos)
                if array[4] == "<DEL>" or array[4] == "<DUP>":
                    pos_end = pos + sv_len
                    self.add_lumpy(chrom, pos_end)
            
            else:
                chrom = array[0]
                pos = int(array[1])

                # if re.search("SVTYPE=BND", line):
                
                fir = re.search("\[(.*?)\[", array[4])
                sec = re.search("](.*?)]", array[4])
                if fir:
                    other_end = fir.group(1)
                else:
                    other_end = sec.group(1)
                locus_info = other_end.split(":")
                chrom2 = locus_info[0]
                pos2 = int (locus_info[1])
                # print (array[4], chrom, pos, abs(pos2 - pos))
                if chrom2 == chrom and abs(pos2 - pos) < min_sv_len:
                    continue
        
                self.add_lumpy(chrom, pos)

    def add_lumpy(self, chrom, pos):
        if chrom in self.all_pos.keys():
            self.all_pos[chrom].append(pos)
        else:
            self.all_pos[chrom]=[pos]

    def cluster_pos(self):
        self.depth_segmentation()
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
        f.close()      

    def get_segments(self):
        for ref in self.ref_len.keys():
            if ref not in self.all_pos:
                self.segments.append([ref, 1, self.ref_len[ref]])
            else:
                start = 1
                for pos in self.all_pos[ref]:
                    if pos - start < MIN_SEG_LEN:
                        start = pos
                        continue
                    self.segments.append([ref, start, pos])
                    start = pos
                if self.ref_len[ref] - start >= MIN_SEG_LEN:
                    self.segments.append([ref, start, self.ref_len[ref]])
        # print (self.segments)
        self.remove_unmapped_segs()
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

    def remove_unmapped_segs_BK(self):
        print ("remove unmapping")
        for i in range(len(self.segments)):
            self.segments[i].append(0)


        for line in open(depth_file):
            array = line.strip().split()
            seg_name = array[0]
            pos = int(array[1])
            dp = int(array[2])
            for i in range(len(self.segments)):
                if self.segments[i][0] == seg_name and pos >= self.segments[i][1] and pos <= self.segments[i][2]:
                    self.segments[i][3] += 1
                    break
        filtered_segments = []
        for i in range(len(self.segments)):
            mapped_ratio = float(self.segments[i][3])/abs(self.segments[i][2]-self.segments[i][1])
            if mapped_ratio > min_mapped_ratio:
                filtered_segments.append(self.segments[i])
            print (self.segments[i], mapped_ratio)
        self.segments = filtered_segments
        print ("unmapping removed")

    def remove_unmapped_segs(self):
        print ("remove unmapping")
        depth_record = {}
        depth_array =[]
        for line in open(depth_file):
            array = line.strip().split()
            seg_name = array[0]
            pos = int(array[1])
            dp = int(array[2])
            if seg_name not in depth_record:
                depth_record[seg_name] = np.zeros(max_chrom_len)
            depth_record[seg_name][pos] = dp
            depth_array.append(dp)
        # min_depth = 0
        # min_depth = np.mean(depth_array) - 5 * np.std(depth_array)
        print ("Median depth:", np.median(depth_array))
        for i in range(len(self.segments)):
            self.segments[i].append(0) #to record hit num
            seg_name = self.segments[i][0]
            start = self.segments[i][1]
            end = self.segments[i][2]
            for pos in range(start, end):
                if seg_name in depth_record and depth_record[seg_name][pos] > 0:
                    self.segments[i][3] += 1
        filtered_segments = []
        for i in range(len(self.segments)):
            mapped_ratio = float(self.segments[i][3])/abs(self.segments[i][2]-self.segments[i][1])
            if mapped_ratio > min_mapped_ratio:
                filtered_segments.append(self.segments[i])
            # print (self.segments[i], mapped_ratio)
        self.segments = filtered_segments
        print ("unmapping removed")

    def depth_segmentation(self):
        intervals = []
        pre_chrom = ""
        start = 0
        index = 0

        f = open(depth_file, 'r')
        for line in f:
            array = line.strip().split()
            chrom = array[0]
            pos = int(array[1])

            if pre_chrom == "":
                start = pos
                index = pos
                pre_chrom = chrom
            elif chrom == pre_chrom:
                if pos - index < min_gap:
                    index = pos
                else:
                    intervals.append([pre_chrom, start, index]) 
                    start = pos
                    index = pos
            else:
                intervals.append([pre_chrom, start, index]) 
                start = pos
                index = pos
                pre_chrom = chrom

        intervals.append([pre_chrom, start, index]) 

        # print (intervals)
        f.close()

        for intes in intervals:
            if abs(intes[2] - intes[1]) < min_exist_len:
                continue
            if intes[0] in self.all_pos:
                self.all_pos[intes[0]] += [intes[1], intes[2]]
            else:
                self.all_pos[intes[0]] = [intes[1], intes[2]]
        # for ref in self.all_pos:
        #     ref_pos_list = self.all_pos[ref]

    def breakpoint_in_bam(self):
        samfile = pysam.AlignmentFile(filename = bam_file, mode = 'rb')
        for read in samfile:
            flag, pos_clip = check_clipped(read)
            if flag:
                chrom = read.reference_name
                pos = read.reference_start + pos_clip
                self.add_lumpy(chrom, pos)


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
lumpy = int(sys.argv[5])
bam_file = sys.argv[6]

bed_file = ref_file + ".bed"


min_gap = 20
MIN_SEG_LEN = 50
min_sv_len = MIN_SEG_LEN
min_exist_len = MIN_SEG_LEN
min_mapped_ratio = 0.8
max_chrom_len = 10000000


print ("seperate reference...")
if lumpy == 0:
    my_bkps = My_bkps()
    my_bkps.read_bkp()
    my_bkps.cluster_pos()
    my_bkps.get_segments()
elif lumpy == 1:
    my_bkps = My_bkps()
    my_bkps.read_lumpy_vcf()
    # my_bkps.breakpoint_in_bam()
    my_bkps.cluster_pos()
    my_bkps.get_segments()
else:
    print ("wrong arg")