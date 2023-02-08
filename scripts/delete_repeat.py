import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re


origin_ref = sys.argv[1]
clean_ref = sys.argv[2]
blast_file = sys.argv[3]
min_repeat_len = 10
end_gap = 20


def read_blast():
    delete_dict = {}
    f = open(blast_file)
    for line in f:
        array = line.strip().split()
        # if array[0] != array[1]:
        if (re.search("NODE_", array[0]) and not re.search("NODE_", array[1])) or\
            (re.search("NODE_", array[1]) and not re.search("NODE_", array[0])):
            chrom = array[1]
            start = int(array[8])
            end = int(array[9])
            if end < start:
                start = int(array[9])
                end = int(array[8])

            if end - start >= min_repeat_len:
                if chrom not in delete_dict:
                    delete_dict[chrom] = [[start, end]]
                else:
                    delete_dict[chrom].append([start, end])
    return delete_dict


def get_clean_interval(delete_interval, seq_len):
    start = 0
    end = seq_len
    for interval in delete_interval:
        if interval[0] < end_gap:
            start = interval[1]
        if seq_len - interval[1] < end_gap:
            end = interval[0]
    return start, end



def main():
    delete_dict = read_blast()
    fasta_sequences = SeqIO.parse(open(origin_ref),'fasta')    
    out_sequences = []
    for record in fasta_sequences:
        if record.id in delete_dict:
            sequence = str(record.seq)
            seq_len = len(sequence)
            delete_interval = delete_dict[record.id]
            start, end = get_clean_interval(delete_interval, seq_len)
            print (record.id, delete_interval, start, end)
            # sequence = sequence[:delete_interval[0]] + sequence[delete_interval[1]:]
            sequence = sequence[start:end]
            record.seq = Seq(sequence)
        out_sequences.append(record)

    SeqIO.write(out_sequences, clean_ref, "fasta")


main()



