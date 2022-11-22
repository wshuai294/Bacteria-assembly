"""
compute N50 based on the blast result
Nov 21, 2022
"""

import sys
from Bio import SeqIO
import re
import csv

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
    # save the mapped interval of each contigs
    # for spliting the original contigs 
    contig_dict = {}
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
        if contig_name not in contig_dict:
            contig_dict[contig_name] = []
        if q_start > q_end:
            a = q_end.copy()
            q_end = q_start
            q_start = a
        contig_dict[contig_name].append([q_start, q_end])
    return contig_dict

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


def get_fasta_len():
    fasta_len = 0
    with open(true) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            # print(record.id)
            fasta_len += len(record.seq)
    return fasta_len

def main():
    list_of_lengths = []
    for contig in contig_dict:
        interval_list = contig_dict[contig]
        print (interval_list)
        interval_list = sort_interval(interval_list)
        print (interval_list)
        interval_list = merge_interval(interval_list)
        print (interval_list)
        for interval in interval_list:
            list_of_lengths.append(interval[1] - interval[0])

    total_match_len = sum(list_of_lengths)
    n50 = calculate_N50(list_of_lengths)
    fasta_len = get_fasta_len()
    completness = round(total_match_len/fasta_len, 2)
    ID = sample.split("/")[-1]
    f = open(f"{sample}.assessment", "w")
    print ("%s\tN50 is %s, Completeness is %s "%(ID, n50, completness))
    print ("%s\tN50 is %s, Completeness is %s "%(ID, n50, completness), file = f)
    f.close()


if __name__ == "__main__":  
    tolerate_gap = 200
    min_map_len = 200
    blast_file = sys.argv[1]
    true = sys.argv[2]
    sample = sys.argv[3]
    contig_dict = read_blast(blast_file)
    main()