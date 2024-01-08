import sys
from Bio import SeqIO
import re
import csv


mummer_alignment = sys.argv[1]
true = sys.argv[2]
sample = sys.argv[3]
new_mummer_alignment = sys.argv[4]


def get_fasta_len():
    fasta_len = 0
    with open(true) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            # print(record.id)
            fasta_len += len(record.seq)
    return fasta_len


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

def read_mummer():
    match_list = []
    f = open(new_mummer_alignment)
    for line in f:
        if line.strip() == '':
            continue
        if line[0] == ">":
            continue
        array = line.strip().split()
        match_len = int(array[-1])
        if match_len < min_match_len:
            continue

        match_list.append(match_len)
    f.close()
    return match_list

def revise_mummer():
    match_list = []
    f = open(mummer_alignment)
    out = open(new_mummer_alignment, 'w')
    writer = csv.writer(out, delimiter ='\t')
    pre_array = []
    
    for line in f:
        line = line.strip()
        if line[0] == ">":
            if len(pre_array) != 0:
                # print (pre_array)
                writer.writerow(pre_array)
                pre_array = []
            array = line.strip().split()
            if array[-1] == "Reverse":
                forward = False
            else:
                forward = True
            print (line, file = out)
            
            continue
        array = line.strip().split()  
        # print (array, pre_array)
        if len(pre_array) != 0:
            if forward:
                gap = int(pre_array[0]) + int(pre_array[-1]) - int(array[0])
            else:
                gap = int(array[0]) - (int(pre_array[0]) - int(pre_array[-1]))
            if gap >= 0:
                pre_array[-1] = int(pre_array[-1])
            elif  gap > tolerate_gap:
                pre_array[-1] = int(pre_array[-1]) + int(array[-1]) 
            else:
                # print (pre_array)
                writer.writerow(pre_array)
                pre_array = array 
        else:
            pre_array = array 
    writer.writerow(pre_array)
    f.close()
    out.close()
    

def assess():
    revise_mummer()
    match_list = read_mummer()
    total_match_len = sum(match_list)
    n50 = calculate_N50(match_list)
    fasta_len = get_fasta_len()
    completness = round(total_match_len/fasta_len, 2)
    ID = sample.split("/")[-1]
    f = open(f"{sample}.assessment", "w")
    print ("%s\tN50 is %s, Completeness is %s "%(ID, n50, completness))
    print ("%s\tN50 is %s, Completeness is %s "%(ID, n50, completness), file = f)
    f.close()

        

if __name__ == "__main__":  
    tolerate_gap = -200
    min_match_len = 1000
    assess()