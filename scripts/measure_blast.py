"""
compute N50 based on the blast and minimap2 results
Nov 21, 2022

PAF formate: https://github.com/lh3/miniasm/blob/master/PAF.md
"""

import sys
from Bio import SeqIO
import re
import csv
import os

def calculate_N50(list_of_lengths):
    """Calculate N50 for a sequence of numbers.
 
    Args:
        list_of_lengths (list): List of numbers.
 
    Returns:
        float: N50 value.
 
    """
    tmp = []
    for tmp_number in set(list_of_lengths):
            tmp += [tmp_number] * list_of_lengths.count(tmp_number) * tmp_number
    tmp.sort()
 
    if (len(tmp) % 2) == 0:
        median = (tmp[int(len(tmp) / 2) - 1] + tmp[int(len(tmp) / 2)]) / 2
    else:
        median = tmp[int(len(tmp) / 2)]
 
    return median

def read_blast(blast_file): # outfmt 7
    blast()
    # save the mapped interval of each contigs
    # for spliting the original contigs 
    contig_dict = {}
    truth_dict = {}
    f = open(blast_file, 'r')
    for line in f:
        if line[0] == "#":
            continue
        array = line.strip().split()
        q_start = int(array[6])
        q_end = int(array[7])
        if abs(q_start - q_end) < min_map_len:
            continue
        contig_name = array[0]
        truth_name = array[1]
        if contig_name not in contig_dict:
            contig_dict[contig_name] = []
        if truth_name not in truth_dict:
            truth_dict[truth_name] = []

        if q_start > q_end:
            a = q_end.copy()
            q_end = q_start
            q_start = a
        contig_dict[contig_name].append([q_start, q_end])

        if int(array[8]) > int(array[9]):
            truth_dict[truth_name].append([int(array[9]), int(array[8])])
        else:
            truth_dict[truth_name].append([int(array[8]), int(array[9])])
    return contig_dict, truth_dict

def sort_interval(interval_list):
    if len(interval_list) < 2:
        return interval_list
    flag = True
    while flag:
        flag = False
        for i in range(1, len(interval_list)):
            if interval_list[i-1][0] > interval_list[i][0]:
                a = interval_list[i-1].copy()
                interval_list[i-1] = interval_list[i]
                interval_list[i] = a
                flag = True
    return interval_list

def merge_interval(interval_list):
    if len(interval_list) < 2:
        return interval_list
    flag = True
    while flag:
        flag = False
        interval_list_copy = interval_list.copy()
        for i in range(1, len(interval_list)):
            if interval_list[i-1][0]-tolerate_gap <= interval_list[i][0] and interval_list[i-1][1] + tolerate_gap >= interval_list[i][0]:
                if interval_list_copy[i][1] > interval_list_copy[i-1][1]:
                    interval_list_copy[i-1][1] = interval_list_copy[i][1]
                del interval_list_copy[i]
                flag = True
                interval_list = interval_list_copy
                break
    return interval_list

def get_fasta_len(fasta):
    fasta_len = 0
    with open(fasta) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            # print(record.id)
            fasta_len += len(record.seq)
    return fasta_len

def main_blast():
    list_of_lengths = []
    for contig in contig_dict:
        interval_list = contig_dict[contig]
        interval_list = sort_interval(interval_list)
        interval_list = merge_interval(interval_list)
        for interval in interval_list:
            list_of_lengths.append(interval[1] - interval[0])

    total_match_len = sum(list_of_lengths)
    n50 = calculate_N50(list_of_lengths)
    result_fasta_len = get_fasta_len(result_fasta)
    precision = round(total_match_len/result_fasta_len, 2)

    true_list_of_lengths = []
    for contig in truth_dict:
        interval_list = truth_dict[contig]
        interval_list = sort_interval(interval_list)
        interval_list = merge_interval(interval_list)
        for interval in interval_list:
            true_list_of_lengths.append(interval[1] - interval[0])
    truth_match_len = sum(true_list_of_lengths)
    truth_n50 = calculate_N50(true_list_of_lengths)
    result_fasta_len = get_fasta_len(result_fasta)
    true_fasta_len = get_fasta_len(true)
    recall = round(truth_match_len/true_fasta_len, 2)
    return ID, round(n50), precision, recall, round(truth_n50)
    


def blast():
    command = """
        true=%s
        result=%s
        sample=%s
        makeblastdb -in $true -dbtype nucl -out $true.db -parse_seqids
        blastn -query $result -db $true.db -outfmt 7 -out $sample.map2true.out
    """%(true, result_fasta, sample)
    os.system(command)

def main_minimap2():
    minimap2_file = sample + ".map2true.paf"
    command = """
    minimap2 --secondary=no %s %s >%s
    """%(true, result_fasta, minimap2_file)
    os.system(command)
    minimap_list_of_lengths = []
    for line in open(minimap2_file):
        array = line.strip().split()
        Alignment_block_length = int(array[9])
        if Alignment_block_length < min_map_len:
            continue
        minimap_list_of_lengths.append(Alignment_block_length)
    print (sorted(minimap_list_of_lengths))
    minimap_n50 = calculate_N50(minimap_list_of_lengths)
    return round(minimap_n50)

def get_direct_N50(result_fasta):
    direct_length_list = []
    with open(result_fasta) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if len(record.seq) < min_map_len:
                continue
            direct_length_list.append(len(record.seq)) 
    direct_N50 = calculate_N50(direct_length_list)
    print (sorted(direct_length_list))
    return round(direct_N50), sum(direct_length_list)
     


if __name__ == "__main__":  
    tolerate_gap = 200
    min_map_len = 200
    true = sys.argv[1]
    sample = sys.argv[2]
    result_fasta = sys.argv[3]
    blast_file = sample + ".map2true.out"
    ID = sample.split("/")[-1]

    contig_dict, truth_dict = read_blast(blast_file)
    ID, n50, precision, recall, truth_n50 = main_blast()
    minimap_n50 = main_minimap2()
    direct_N50, result_size = get_direct_N50(result_fasta) 
    f = open(f"{sample}.assessment", "w")
    print ("%s\tPrecision: %s; Recall: %s; Origin N50: %s, Minimap_N50: %s, size:%s"%(ID, precision, recall, direct_N50, minimap_n50,result_size))
    print ("%s\tPrecision: %s; Recall: %s;Origin N50: %s, Minimap_N50: %s,size:%s;"%(ID, precision, recall, direct_N50, minimap_n50,result_size), file = f)
    f.close()

    