import sys
from Bio import SeqIO


mummer_alignment = sys.argv[1]
true = sys.argv[2]
sample = sys.argv[3]


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
    f = open(mummer_alignment)
    for line in f:
        if line[0] == ">":
            continue
        array = line.strip().split()
        match_len = int(array[-1])
        if match_len < min_match_len:
            continue

        match_list.append(match_len)
    f.close()
    return match_list

def assess():
    match_list = read_mummer()
    total_match_len = sum(match_list)
    n50 = calculate_N50(match_list)
    fasta_len = get_fasta_len()
    completness = round(total_match_len/fasta_len, 2)
    f = open(f"{sample}.assessment", "w")
    print ("N50 is %s, Completeness is %s "%(n50, completness))
    print ("N50 is %s, Completeness is %s "%(n50, completness), file = f)
    f.close()

        

    

min_match_len = 10000
assess()