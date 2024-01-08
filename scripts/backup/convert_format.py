import sys
from Bio import SeqIO
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re



original_graph = sys.argv[1]
out_path = sys.argv[2]
seg_ref = sys.argv[3]

# original_graph = "/home/wangshuai/assembly_result/real_data/DRR198804.graph.txt"
# out_path = "/home/wangshuai/assembly_result/real_data/"
# seg_ref = "/home/wangshuai/assembly_result/real_data/DRR198804.seg.fa"



def get_contig_lengths(fasta_file):
    contig_lengths = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        contig_lengths[record.id] = len(record.seq)
    return contig_lengths


contig_lengths = get_contig_lengths(seg_ref)
seg_out = out_path + "/input.segs"
junc_out = out_path + "/input.juncs"


out = open(seg_out, "w")
jout = open(junc_out, "w")
seg_index = 0
record_index = {}
for line in open(original_graph):
    array = line.strip().split()
    if array[0] == "SEG":
        seg_index += 1
        record_index[array[1]] = seg_index

        
        print (seg_index, array[1], 1, contig_lengths[array[1]], sep = "\t", file = out)
    
    elif array[0] == "JUNC":
        print (record_index[array[1]], array[2], record_index[array[3]], array[4], sep = "\t", file = jout)
        
out.close()
jout.close()

