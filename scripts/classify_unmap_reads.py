"""
unmapped reads are mixture of paired and single reads
classify them
"""

import sys


raw_fq1 = sys.argv[1]
raw_fq2 = sys.argv[2]
outdir = sys.argv[3]
ID = sys.argv[4]
fq1 = outdir + "/" + ID + ".unmapped.1.fq"
fq2 = outdir + "/" + ID + ".unmapped.2.fq"
fqs = outdir + "/" + ID + ".unmapped.s.fq"


def extract_read_name(fastq):
    read_name_dict = {}
    i = 0
    f = open(fastq, 'r') 
    for line in f:
        if i % 4 == 0:
            read_name = line[1:].strip().split("/")[0]
            read_name_dict[read_name] = 1
        i += 1
    f.close()
    return read_name_dict

def extract_single_reads(read_name_dict, fastq, clean_fastq, single_fastq, flag):
    clean_f = open(clean_fastq, 'w')
    if flag == "initial":
        single_f = open(single_fastq, 'w')
    else:
        single_f = open(single_fastq, 'w+')
    if_paired = False
    f = open(fastq, 'r') 
    i = 0
    for line in f:
        if i % 4 == 0:
            read_name = line[1:].strip().split("/")[0]
            if read_name in read_name_dict:
                # print ("*********************")
                if_paired = True
            else:
                if_paired = False
        if if_paired:
            print (line, end = '', file = clean_f)
        else:
            print (line, end = '', file = single_f)
        i += 1
    clean_f.close()
    single_f.close()

def main():
    fq1_read_name_dict = extract_read_name(raw_fq1)
    extract_single_reads(fq1_read_name_dict, raw_fq2, fq2, fqs, "initial")
    # print (fq1_read_name_dict)
    fq2_read_name_dict = extract_read_name(raw_fq2)
    extract_single_reads(fq2_read_name_dict, raw_fq1, fq1, fqs, "not initial")

main()