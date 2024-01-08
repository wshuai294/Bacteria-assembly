from Bio import SeqIO
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re


original_ref = sys.argv[1]
graph_file = sys.argv[2]
contigs = sys.argv[3]


def get_contig_lengths(fasta_file):
    contig_lengths = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        contig_lengths[record.id] = len(record.seq)
    return contig_lengths


def main():
    used_segs = {}
    output_handle = open(contigs, "w")
    record_dict = SeqIO.to_dict(SeqIO.parse(original_ref, "fasta"))
    f = open(graph_file, 'r')
    contig_index = 1
    for line in f:
        if line.strip() == '':
            continue
        segs_list = line.strip().split()

        # if segs_list[0][-1] == "'": # check if the path is already recorded.
        #     seg_name = segs_list[0][:-1]
        # else:
        #     seg_name = segs_list[0]
        # if seg_name in used_segs:
        #     continue

        for i in range(len(segs_list)):
            seg_name = segs_list[i][:-1]
            record = record_dict[seg_name]
            if segs_list[i][-1] == "-":
                # add_seq = str(record.seq)[::-1]
                add_seq = str(record.seq.reverse_complement())
            else:
                add_seq = str(record.seq)
                # segs_list[i] = segs_list[i][:-1]

                # record.seq = record.seq.reverse_complement()
            # else:
            #     seg_name = segs_list[i][:-1]
            #     record = record_dict[seg_name]
            used_segs[seg_name] = 1
            if i == 0:
                merged_record = SeqRecord(Seq(""), f"contig_{contig_index}", '', '')
            merged_record.seq = Seq(str(merged_record.seq) + add_seq)
            merged_record.description += "\t" + segs_list[i]
        SeqIO.write(merged_record, output_handle, "fasta")
        # print (f"contig_{contig_index}", len(merged_record.seq))   
        contig_index += 1

    output_handle.close()
    f.close()

def compute_node_length():
    f = open(graph_file, 'r')
    h = open(graph_file + ".length", 'w')
    y = open(graph_file + ".ref.bed", 'w')
    contig_lengths = get_contig_lengths(original_ref)

    contig_index = 1
    for line in f:
        line = line.strip()
        if line == '':
            continue
        
            print (line, file = h)

        segs_list = line.strip().split()
        accumu_len = 0
        for i in range(len(segs_list)):
            seg_name = segs_list[i][:-1]
            if re.search("NODE_", seg_name):
                # length = int(seg_name.split("_")[3])
                length = contig_lengths[seg_name]
                new_seg_name = seg_name
            else:
                array = seg_name.split(":")[-1].split("-")
                # length = abs(int(array[1]) - int(array[0]))
                length = contig_lengths[seg_name]
                new_seg_name = seg_name + "_len_" + str(length)
                print ("contig_%s\t%s\t%s"%(contig_index, accumu_len, accumu_len+length), file = y)
            accumu_len += length
            new_seg_name = new_seg_name + "_accu_" + str(accumu_len)
            print (new_seg_name, end = "\t", file = h)
        print ('', end = "\n", file = h)

        contig_index += 1

    f.close()
    h.close()
    y.close()


main()
compute_node_length()