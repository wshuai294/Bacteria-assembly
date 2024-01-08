import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re


origin_ref = sys.argv[1]
clean_ref = sys.argv[2]
blast_file = sys.argv[3]
min_repeat_len = 100
end_gap = 20
repeat_cutoff = 0.95


def read_blast():
    delete_dict = {}
    delete_segment = {}
    kept_segment = {}
    
    f = open(blast_file)
    for line in f:
        array = line.strip().split()
        identity = float(array[2])
        if identity < 95:
            continue

        chrom = array[1]
        start = int(array[8])
        end = int(array[9])
        if end < start:
            start = int(array[9])
            end = int(array[8])

        if array[0] == array[1]:
            continue
        elif re.search("NODE_", array[1]) and not re.search("NODE_", array[0]) and not re.search("NODE_INS", array[1]): #  (re.search("NODE_", array[0]) and not re.search("NODE_", array[1]))

            if end - start >= min_repeat_len:
                # print (array)
                if chrom not in delete_dict:
                    delete_dict[chrom] = [[start, end]]
                else:
                    delete_dict[chrom].append([start, end])
        ### keep one representative seg for repeat segments
        elif not re.search("NODE_", array[0]) and not re.search("NODE_", array[1]):
            seg_1_len = cal_seg_len(array[0])
            seg_2_len = cal_seg_len(array[1])
            match_len = int(array[3])
            if match_len/seg_1_len > repeat_cutoff and match_len/seg_2_len > repeat_cutoff:
                print ("%s and %s are repeats."%(array[0], array[1]))
                # if array[1] not in kept_segment:
                delete_segment[array[1]] = 1
                kept_segment[array[0]] = 1
                
            else:
                if end - start >= min_repeat_len:
                    # print (array)
                    if chrom not in delete_dict:
                        delete_dict[chrom] = [[start, end]]
                    else:
                        delete_dict[chrom].append([start, end])

    # delete_segment = {"NZ_CP043539.1:717052-729238": 1, "NZ_CP043539.1:5180261-5191067": 1}
    return delete_dict, delete_segment

def cal_seg_len(seg):
    # example NZ_CP043539.1:717052-729238
    new_array = seg.split(":")[1].split("-")
    # print (new_array)
    length = abs(int(new_array[0]) - int(new_array[1]))  
    return length

def get_clean_interval_bk(delete_interval, seq_len):
    start = 0
    end = seq_len
    for interval in delete_interval:
        if interval[0] < end_gap:
            start = interval[1]
        if seq_len - interval[1] < end_gap:
            end = interval[0]
    if end > start:
        return start, end, "keep"
    else:
        return start, end, "remove"

def merge_intervals(intervals):
    sorted_intervals = sorted(intervals, key=lambda x: x[0])
    merged_intervals = []
    for interval in sorted_intervals:
        if not merged_intervals or interval[0] > merged_intervals[-1][1]:
            merged_intervals.append(interval)
        else:
            merged_intervals[-1] = (merged_intervals[-1][0], max(merged_intervals[-1][1], interval[1]))
    return merged_intervals

def get_clean_interval(delete_interval, seq_len):
    merged_intervals = merge_intervals(delete_interval)
    clean_intervals = []

    # iterate over the given list of intervals, removing each interval from the list of all possible intervals
    start = 0
    for interval in merged_intervals:
        end = interval[0]
        if end - start > 1000:
            clean_intervals.append([start, end])
        start = interval[1]
    end = seq_len
    if end - start > 1000:
        clean_intervals.append([start, end])

    # print ("repeat interval", merged_intervals)
    # print ("clean interval", clean_intervals)
    if len(clean_intervals) > 0:
        flag = "keep"
    else:
        flag = "remove"
    return clean_intervals, flag

def sort_and_merge(intervals):
    # Sort the intervals by their start points
    intervals = sorted(intervals, key=lambda x: x[0])
    
    # Initialize the current interval and result list
    current_interval = intervals[0]
    result = []
    
    # Iterate over the remaining intervals
    for interval in intervals[1:]:
        # If the current interval overlaps with the next interval, merge them
        if current_interval[1] >= interval[0]:
            current_interval = (current_interval[0], max(current_interval[1], interval[1]))
        # Otherwise, add the current interval to the result list and set the current interval to the next interval
        else:
            result.append(current_interval)
            current_interval = interval
    
    # Add the final current interval to the result list
    result.append(current_interval)
    delete_length = get_total_length(result)
    return result, delete_length

def get_total_length(intervals):
    total_length = 0
    for interval in intervals:
        total_length += interval[1] - interval[0]
    return total_length

def main():
    delete_dict, delete_segment = read_blast()
    # print (delete_dict)
    print ("<<<<< delete_segment num", len(delete_segment))
    fasta_sequences = SeqIO.parse(open(origin_ref),'fasta')    
    out_sequences = []
    total_delete_length = 0
    total_remain_length = 0
    remain_contig_num = 0
    for record in fasta_sequences:

        if record.id in delete_segment:
            seg_len = cal_seg_len(record.id)
            total_delete_length += seg_len
            print ("4. remove this seg", record.id, "len is", seg_len)
            continue
            
        elif record.id in delete_dict:
            sequence = str(record.seq)
            seq_len = len(sequence)
            delete_interval = delete_dict[record.id]
            merged_delete_interval, delete_length = sort_and_merge(delete_interval)
            total_delete_length += delete_length
            if seq_len - delete_length < 20:
                print ("1. remove this seg", record.id, seq_len)
                continue
            clean_intervals, flag = get_clean_interval(delete_interval, seq_len)
            if flag == "remove":
                print ("2. remove this seg", record.id, seq_len)
                continue
            print ("3. keep parts", record.id)
            total_remain_length += get_total_length(clean_intervals)
            for clean in clean_intervals:
                # sequence = sequence[:delete_interval[0]] + sequence[delete_interval[1]:]
                name = record.id + "%" + str(clean[0]) + "%" + str(clean[1])
                part_sequence = sequence[clean[0]:clean[1]]
                # print (name, sequence)
                part_record = SeqRecord(Seq(part_sequence), id=name, description="")
                # record.seq = Seq(sequence)
                out_sequences.append(part_record)
                remain_contig_num += 1
            continue
        else:
            total_remain_length += len(record.seq)
            remain_contig_num += 1

        out_sequences.append(record)
    print ("final delete length:", total_delete_length, "remain length:", total_remain_length, "remain contig:", remain_contig_num)
    SeqIO.write(out_sequences, clean_ref, "fasta")


main()



