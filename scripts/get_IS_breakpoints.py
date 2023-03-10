"""
map reference assembly to IS database
find the IS breakpints
"""

import sys, os

# ref_assembly = sys.argv[1]

IS_DB = "/home/wangshuai/assembly_result/IS.fna" 
### https://github.com/thanhleviet/ISfinder-sequences/blob/master/IS.fna

def blast(ref_assembly, blast_file):
    command = f"blastn -query {ref_assembly} -subject {IS_DB} -out {blast_file} -outfmt 6"
    os.system(command)

    

def get_bkp(blast_file):
    bkp_list = []
    f = open(blast_file)
    for line in f:
        array = line.strip().split()
        identity = float(array[2])
        match_len = int(array[3])
        chrom = array[0]
        start = int(array[6])
        end = int(array[7])
        if identity > 99 and match_len > 500:
            bkp_list.append([chrom, start])
            bkp_list.append([chrom, end])
    f.close()
    # print (bkp_list)
    return bkp_list

def IS_bkp(ref_assembly, sample):
    blast_file = f"{sample}.blast.IS.result"
    blast(ref_assembly, blast_file)
    bkp_list = get_bkp(blast_file)
    return bkp_list


def get_repeat_bkp(ref_assembly, sample):
    blast_file = f"{sample}.blast.self.result"
    command = f"blastn -query {ref_assembly} -subject {ref_assembly} -out {blast_file} -outfmt 6"
    os.system(command)
    bkp_list = []
    f = open(blast_file)
    for line in f:
        array = line.strip().split()
        identity = float(array[2])
        match_len = int(array[3])
        chrom = array[0]
        start = int(array[6])
        end = int(array[7])
        if identity > 99 and match_len > 2000:
            bkp_list.append([chrom, start])
            bkp_list.append([chrom, end])
    f.close()
    # print (bkp_list)
    return bkp_list


# test = "/home/wangshuai/assembly_result/test.out"
# get_bkp(test)


