from Bio import SeqIO
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


original_ref = sys.argv[1]
graph_file = sys.argv[2]
contigs = sys.argv[3]


def main():
    used_segs = {}
    output_handle = open(contigs, "w")
    record_dict = SeqIO.to_dict(SeqIO.parse(original_ref, "fasta"))
    f = open(graph_file, 'r')
    contig_index = 1
    for line in f:
        segs_list = line.strip().split()

        # if segs_list[0][-1] == "'": # check if the path is already recorded.
        #     seg_name = segs_list[0][:-1]
        # else:
        #     seg_name = segs_list[0]
        # if seg_name in used_segs:
        #     continue

        for i in range(len(segs_list)):
            if segs_list[i][-1] == "'":
                # segs_list[i] = segs_list[i][:-1]
                seg_name = segs_list[i][:-1]
                record = record_dict[segs_list[i][:-1]]
                record.seq = record.seq.reverse_complement()
            else:
                seg_name = segs_list[i]
                record = record_dict[segs_list[i]]
            used_segs[seg_name] = 1
            if i == 0:
                merged_record = SeqRecord(Seq(""), f"contig_{contig_index}", '', '')

            merged_record.seq = Seq(str(merged_record.seq) + str(record.seq))
            merged_record.description += "\t" + segs_list[i]
        SeqIO.write(merged_record, output_handle, "fasta")
        print (f"contig_{contig_index}", len(merged_record.seq))   
        contig_index += 1

    output_handle.close()
    f.close()


main()