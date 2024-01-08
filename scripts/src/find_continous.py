from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import sys

def split_fasta(fasta_file, interval_list, output_file):
    output_records = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id
        sequence = str(record.seq)

        for i, (start, end) in enumerate(interval_list, start=1):
            if end - start < minimum_seg_len:
                continue
            split_sequence = sequence[start-1:end]
            split_id = f"{seq_id}:{start}-{end}"
            split_record = SeqRecord(Seq(split_sequence), id=split_id, description="")
            output_records.append(split_record)

    SeqIO.write(output_records, output_file, "fasta")


def merge_intervals(intervals):
    min_overlap_len = 1
    intervals.sort(key=lambda x: x[0])  # Sort intervals based on start time
    merged = [intervals[0]]

    for interval in intervals[1:]:
        # if interval[0] <= merged[-1][1]:
        if merged[-1][1] - interval[0] > min_overlap_len:  # Overlapping intervals
            merged[-1] = (merged[-1][0], max(merged[-1][1], interval[1]))
        else:  # Non-overlapping interval
            merged.append(interval)

    return merged


def collect_intervals(file):
    map_intervals = []

    my_read = -1
    my_interval = [float("inf"), 0]
    my_read_interval = [-1, -1]
    start_array = [-1, -1, False, -1]
    cover_flag = False
    cover_interval = [-1, -1]

    for line in open(file):
        array = line.strip().split()
        line_index = int(array[0])
        read_pos = int(array[1])
        coder_index = int(array[2])
        ref_pos = int(array[3])
        kmer_ambiguity = int(array[4])

        # if kmer_ambiguity > 1:
        #     continue
        if ref_pos == 0:
            continue
        if line_index != my_read:
            if my_read_interval[0] != -1:

                # read_diff = abs(my_read_interval[1] - my_read_interval[0])
                # ref_diff = abs(my_interval[1] - my_interval[0])

                # if abs(read_diff - ref_diff) < 30:
                #     map_intervals.append(my_interval) 
                #     test_flag = True
                # else:
                #     test_flag = False
                # #     
                #     print ("##yes ", my_read, my_interval, my_read_interval)
                # else:
                #     print (my_read, my_interval, my_read_interval)
                if cover_flag:
                    map_intervals.append(cover_interval) 
                # if cover_flag and not test_flag:
                #     map_intervals.append(cover_interval) 
                #     # if cover_interval[0] != my_interval[0] or cover_interval[1] != my_interval[1]:
                #     print (my_read, my_interval, my_read_interval, cover_flag, cover_interval, test_flag)
                    # print ("final", cover_interval)
                # else:
                #     print (my_read, my_interval, my_read_interval, cover_flag, cover_interval)

            my_read = line_index
            my_interval = [-1, -1]
            my_read_interval = [-1, -1]
            cover_flag = False
            cover_interval = [-1, -1]

        
        # if read_pos > 30:
        if True:
            read_diff = abs(my_read_interval[1] - my_read_interval[0])
            ref_diff = abs(my_interval[1] - my_interval[0])

            if abs(read_diff - ref_diff) < 30:# and ref_diff > 30:
                cover_flag = True
                cover_interval = my_interval.copy()

        # if my_read == 365293:
        #     print (my_read, my_interval, my_read_interval, cover_flag, cover_interval, read_diff, ref_diff, abs(read_diff - ref_diff))
        
        if coder_index == 0:
            start_array = [read_pos, ref_pos, True, 0]
        else:
            if ref_pos != start_array[1]:  # three coders should have the same ref pos
                start_array[2] = False
        if coder_index == 2 and start_array[2]:
            # print (my_read, "start", start_array)
            # break
            if my_interval[0] == -1:
                my_interval = [start_array[1], start_array[1]]
            if start_array[1] <= my_interval[0]:
                my_interval[0] = start_array[1]
                my_read_interval[0] = start_array[0]
            if start_array[1] >= my_interval[1]:
                my_interval[1] = start_array[1]   
                my_read_interval[1] = start_array[0]   
    return map_intervals


if __name__ == "__main__":
    # ref_file = "/mnt/d/breakpoints/assembly/sim/database/Escherichia_coli/NC_007779.1.fasta"
    # file = "/home/wangshuai/assembly_result/test.map.tab"
    # out_file = "/home/wangshuai/assembly_result/test.split.fasta"

    ref_file = sys.argv[1]
    file = sys.argv[2]
    out_file = sys.argv[3]
    minimum_seg_len = int(sys.argv[4])


    # ref_file = "/mnt/d/breakpoints/assembly/sim/database/Escherichia_coli/NZ_CP028685.1.fasta"  #803  804
    # truth = "/mnt/d/breakpoints/assembly/sim/standard/truth/GCF_000008865.2_ASM886v2_genomic.fna"
    # truth = "/mnt/d/breakpoints/assembly/sim/standard/truth/GCF_000005845.2_ASM584v2_genomic.fna"  #805  806



    print ("find continous intervals...")
    map_intervals = collect_intervals(file)
    print ("collected interval count", len(map_intervals))
    merged = merge_intervals(map_intervals)
    total_len = 0
    for interval in merged:
        total_len += (interval[1] - interval[0])
    # print (map_intervals)
    # print (merged)
    split_fasta(ref_file, merged, out_file)
    print (total_len, len(merged))


    # command = f""" rm -r /home/wangshuai/assembly_result/test
    # ~/softwares/quast-5.2.0/quast.py -t 10 {out_file} -r {truth} -o /home/wangshuai/assembly_result/test --strict-NA 
    #           cat /home/wangshuai/assembly_result/test/report.txt"""

    # os.system(command)

    