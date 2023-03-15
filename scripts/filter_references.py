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
    IS_dict = {}
    f = open(blast_file)
    for line in f:
        array = line.strip().split()
        identity = float(array[2])
        match_len = int(array[3])
        chrom = array[0]
        start = int(array[6])
        end = int(array[7])
        IS = array[1]
        if identity > 99 and match_len > 700:
            if IS not in IS_dict:
                IS_dict[IS] = 0
            IS_dict[IS] += 1
            # bkp_list.append([chrom, start])
            # bkp_list.append([chrom, end])
    f.close()
    sorted_IS_dict = sorted(IS_dict.items(), key=lambda item: item[1], reverse = True)
    # print (bkp_list)
    return sorted_IS_dict

# def IS_bkp(ref_assembly):
#     blast_file = f"{sample}.blast.IS.result"
#     blast(ref_assembly, blast_file)
#     bkp_list = get_bkp(blast_file)
#     return bkp_list


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


def count_IS_in_each_ref():
    assembly_dict = {}
    ref_list = "/home/wangshuai/assembly_result/self_fasta.txt"
    f = open(ref_list, 'r')
    for line in f:
        ref_assembly = line.strip()
        if not os.path.isfile(ref_assembly):
            continue
        blast_file = ref_assembly + ".IS.blast"
        blast(ref_assembly, blast_file)
        if not os.path.isfile(blast_file):
            continue
        sorted_IS_dict = get_bkp(blast_file)
        # if ref_assembly == "/mnt/d/breakpoints/assembly/sim/database/ecoli/104.fasta":
            # print (ref_assembly, len(sorted_IS_dict), sorted_IS_dict[0][1])
        assembly_dict[ref_assembly] = sorted_IS_dict[0][1]
        # break
    sorted_assembly_dict = sorted(assembly_dict.items(), key=lambda item: item[1])
    # print (sorted_assembly_dict[:100])
    for i in range(500):
        print (sorted_assembly_dict[i][0])

# test = "/home/wangshuai/assembly_result/test.out"
# get_bkp(test)
count_IS_in_each_ref()


