"""
Given reference and sequencing reads, find the reference segments that might present in the sample using kmer. 

"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import sys
import argparse
import numpy as np

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

def collect_intervals_bk(file): # consider kmer map and extract read map region
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

        if kmer_ambiguity > 1:  # important
            continue
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

def collect_intervals(file): ## direct read map interval
    map_intervals = []
    for line in open(file):
        array = line.strip().split()
        array[0] = int(array[0])
        array[1] = int(array[1])
        map_intervals.append(array)

    return map_intervals


def kmer_process():
    command  = f"""
    rm {kmer_count_file}_*.part.txt
    {sys.path[0]}/src/continuous {ref_file} {options.fq1} {options.fq2} {kmer_count_file} {best_sample_ratio} {options.k} {options.t}
    cat {kmer_count_file}_*.part.txt > {kmer_count_file}
    rm {kmer_count_file}_*.part.txt
    """
    # print (command)
    os.system(command)

def check_input():
    if not os.path.exists(options.o):
        os.system(f"mkdir {options.o}")
    if not os.path.exists(options.r):
        print ("Reference file: %s is not detected."%(options.r))
        sys.exit()
    if not os.path.exists(options.fq1):
        print ("Reference file: %s is not detected."%(options.fq1))
        sys.exit()

def count_bases_fasta(fasta_file):
    total_bases = 0

    with open(fasta_file, "r") as fasta_handle:
        for record in SeqIO.parse(fasta_handle, "fasta"):
            total_bases += len(record.seq)

    return total_bases

def count_bases_fastq(fastq_file):
    count = 0
    read_lens = []
    with open(fastq_file, 'r') as file:
        for line_num, line in enumerate(file):
            if line_num % 4 == 1:  # Check the second line of each record
                count += len(line.strip())
                read_lens.append(len(line.strip()))
    print ("read length mean is %s, median is %s."%(np.mean(read_lens), np.median(read_lens)))
    return count

def estimate_sample_ratio(fasta_file, fastq_file):
    # best depth 16x
    est_depth = float(count_bases_fastq(fastq_file))/count_bases_fasta(fasta_file)
    best_sample_ratio = round(options.sample_depth/est_depth * 100)
    return best_sample_ratio



if __name__ == "__main__":


    parser = argparse.ArgumentParser(description="Detect HGT breakpoints from metagenomics sequencing data.", add_help=False, \
    usage="%(prog)s -h", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")
    required.add_argument("-r", type=str, help="<str> best-match reference file.", metavar="\b")
    required.add_argument("--fq1", type=str, help="<str> unzipped fastq 1 file.", metavar="\b")
    required.add_argument("--fq2", type=str, help="<str> unzipped fastq 2 file.", metavar="\b")
    required.add_argument("-s", type=str, default="sample", help="<str> Sample name.", metavar="\b")
    required.add_argument("-o", type=str, default="./", help="<str> Output folder.", metavar="\b")
    optional.add_argument("-m", type=int, default=1000, help="<int> minimum segment length.", metavar="\b")
    optional.add_argument("-k", type=int, default=26, help="<int> kmer length.", metavar="\b")
    optional.add_argument("-t", type=int, default=10, help="<int> number of threads.", metavar="\b")
    # optional.add_argument("--sample_depth", type=float, default=16, help="<float> only retain reads with this depth in downsampling.", metavar="\b")
    # optional.add_argument("--sample_ratio", type=int, default=100, help="<0-100> sample ratio (%).", metavar="\b")
    optional.add_argument("-h", "--help", action="help")

    options = parser.parse_args()

    if len(sys.argv) == 1:
        print (f"see python {sys.argv[0]} -h")
    else:
        check_input()
        ref_file = options.r
        # file = sys.argv[2]
        # out_file = sys.argv[3]
        # minimum_seg_len = int(sys.argv[4])

        kmer_count_file = options.o + "/" + options.s + ".kmer.txt"
        minimum_seg_len = options.m
        out_file = options.o + "/" + options.s + ".split.fasta"
        # best_sample_ratio = estimate_sample_ratio(options.r, options.fq1)
        # print (f"sampling ratio is {best_sample_ratio}.")
        best_sample_ratio = 101 # use all reads

        print ("start kmer matching...")
        kmer_process()
        print ("find reference segments...")
        map_intervals = collect_intervals(kmer_count_file)
        print ("collected interval count", len(map_intervals))
        merged = merge_intervals(map_intervals)
        total_len = 0
        for interval in merged:
            total_len += (interval[1] - interval[0])

        split_fasta(ref_file, merged, out_file)
        # print (total_len, len(merged))
        print ("extract %s segments, total length is %s bp."%(len(merged), total_len))