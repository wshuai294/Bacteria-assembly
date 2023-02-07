"""
unmapped reads are mixture of paired and single reads
classify them
"""

import sys
import os


raw_fq1 = sys.argv[1]
raw_fq2 = sys.argv[2]
unmap_reads_file_fq1 = f"{raw_fq1}.unmap.reads.txt"
unmap_reads_file_fq2 = f"{raw_fq2}.unmap.reads.txt"
outdir = sys.argv[3]
ID = sys.argv[4]
fq1 = outdir + "/" + ID + ".unmapped.1.fq"
fq2 = outdir + "/" + ID + ".unmapped.2.fq"
fqs = outdir + "/" + ID + ".unmapped.s.fq"

def combine_frag():
    # the c++ script genetates _part_ files for each thread
    # combine them and delete temp files.
    cmd = f"""cat {raw_fq1}_part_*.txt >{unmap_reads_file_fq1}
            rm {raw_fq1}_part_*.txt"""
    os.system(cmd)

    cmd = f"""cat {raw_fq2}_part_*.txt >{unmap_reads_file_fq2}
            rm {raw_fq2}_part_*.txt"""
    os.system(cmd)


def extract_read_name(unmap_reads_file):
    read_name_dict = {}
    f = open(unmap_reads_file, 'r') 
    for line in f:
        read_name = line.strip()
        read_name_dict[read_name] = 1
    f.close()
    return read_name_dict

def extract_single_reads(focus_read_name_dict, other_read_name_dict, raw_fastq, clean_fastq, single_fastq, flag):
    clean_f = open(clean_fastq, 'w')
    if flag == "initial":
        single_f = open(single_fastq, 'w')
    else:
        single_f = open(single_fastq, 'w+')
    if_paired = False
    f = open(raw_fastq, 'r') 
    i = 0
    for line in f:
        if i % 4 == 0:
            read_name = line[1:].strip().split("/")[0]
             
            if read_name in other_read_name_dict:
                # 
                if_paired = True
            else:
                if_paired = False
        if read_name in focus_read_name_dict:
            if if_paired:
                print (line, end = '', file = clean_f)
            else:
                print (line, end = '', file = single_f)
        i += 1
    clean_f.close()
    single_f.close()
    f.close()

def main():
    combine_frag()
    fq1_read_name_dict = extract_read_name(unmap_reads_file_fq1)
    # print (fq1_read_name_dict)
    fq2_read_name_dict = extract_read_name(unmap_reads_file_fq2)
    extract_single_reads(fq1_read_name_dict, fq2_read_name_dict, raw_fq1, fq1, fqs, "initial")
    extract_single_reads(fq2_read_name_dict, fq1_read_name_dict, raw_fq2, fq2, fqs, "not initial")

main()