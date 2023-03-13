from Bio import SeqIO
import pysam
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os, sys
import argparse



def jellyfish():
    command = f"""
        jellyfish count -m {k} -s 2000000000 -t {threads} -C {fq1} -o {jf_file}
        jellyfish query {jf_file} -s {input_fasta} > {count_file}
    """
    os.system(command)

def get_contig_interval(fai):
    # if not os.path.isfile(fai):
    os.system(f"samtools faidx {input_fasta}")
    contig_interval = []
    f = open(fai, 'r')
    start = 1
    for line in f:
        array = line.strip().split()   
        contig_name = array[0]
        contig_length = int(array[1])
        end = start + contig_length
        contig_interval.append([start, end-1, contig_name])
        start = end
    f.close()
    return contig_interval

def locus_transfer(intervals, value):
    # Initialize pointers
    left = 0
    right = len(intervals) - 1
    
    # Binary search for interval containing value
    while left <= right:
        mid = (left + right) // 2
        if value >= intervals[mid][0] and value <= intervals[mid][1]:
            contig = intervals[mid][2]
            relative_locus = value - intervals[mid][0] + 1
            return contig, relative_locus
        elif value < intervals[mid][0]:
            right = mid - 1
        else:
            left = mid + 1

def read_kmer(contig_interval):
    unmapped_loci = {}
    index = 1
    f = open(count_file, 'r')
    for line in f:
        array = line.strip().split()
        count = int(array[1])
        if count < min_count:
            # print (index, contig_interval)
            real_locus = index + k -1
            contig, relative_locus = locus_transfer(contig_interval, real_locus)
            if contig not in unmapped_loci:
                unmapped_loci[contig] = []
            unmapped_loci[contig].append(relative_locus)
            # print (index, contig, relative_locus, count)

        index += 1
    f.close()
    return unmapped_loci

def get_breakpoints(unmapped_loci):
    my_breakpoints = {}
    for contig in unmapped_loci:
        loci = unmapped_loci[contig]
        breakpoints = []
        for locus in loci:
            if len(breakpoints) == 0:
                breakpoints.append(locus)
            elif locus - breakpoints[-1] > k:
                breakpoints.append(locus)
        my_breakpoints[contig] = breakpoints
    return my_breakpoints

def split_fasta(input_fasta, output_fasta, my_breakpoints):
    # Load the FASTA file
    fasta_file = SeqIO.parse(input_fasta, "fasta")
    handle = open(output_fasta, "w")
    # Split the FASTA file based on the breakpoints
    for record in fasta_file:
        start = 0
        end = len(record)
        if record.id in my_breakpoints:
            breakpoint_list = my_breakpoints[record.id]
        else:
            breakpoint_list = []
        for breakpoint in breakpoint_list:
            if abs(start - breakpoint) > min_contig_len:
                new_record = record[start:breakpoint]
                new_record.id = f"{record.id}_{start}_{breakpoint}"
                new_record.description = ""
                SeqIO.write(new_record, handle, "fasta")
                start = breakpoint
        if abs(start - end) > min_contig_len:
            new_record = record[start:end]
            new_record.id = f"{record.id}_{start}_{end}"
            new_record.description = ""
            SeqIO.write(new_record, handle, "fasta")


if __name__ == "__main__":
    # k = 14
    # min_count = 10
    # min_contig_len = 200
    # threads = 10

    # input_fasta = "/home/wangshuai/assembly_result/test_short/test_short_ID.contigs.polish_1.fasta"
    # output_fasta = "/home/wangshuai/assembly_result/test_short/test_short_ID.contigs.polish_2.fasta"
    # fq1 = "/home/wangshuai/assembly_result/test_short.1.fq"


    parser = argparse.ArgumentParser(description="Polish assembly with kmer", add_help=False, \
    usage="python3 %(prog)s -h", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("Required arguments")
    optional = parser.add_argument_group("Optional arguments")
    required.add_argument("-i", type=str, help="Input assembly (fasta)", metavar="\b")
    required.add_argument("-r", type=str, help="The fastq1 file", metavar="\b")
    required.add_argument("-o", type=str, help="The polished assembly (fasta).", metavar="\b")
    optional.add_argument("-c", type=int, help="Minimum kmer count", metavar="\b", default=10)
    optional.add_argument("-k", type=int, help="Kmer length", metavar="\b", default=14)
    optional.add_argument("-j", type=int, help="Number of threads.", metavar="\b", default=10)
    optional.add_argument("-e", type=int, help="Minimum contig length", metavar="\b", default=200)
    optional.add_argument("-h", "--help", action="help")
    args = vars(parser.parse_args()) 

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)

    input_fasta = args["i"]
    output_fasta = args["o"]
    fq1 = args["r"]

    k = args["k"]
    min_count = args["c"]
    min_contig_len = args["e"]
    threads = args["j"]




    fai = input_fasta + ".fai"
    jf_file = input_fasta + ".jf"
    count_file = input_fasta + ".count"


    jellyfish()
    contig_interval = get_contig_interval(fai)
    # print (contig_interval)
    # contig, relative_locus = locus_transfer(contig_interval, 270931)
    # print (contig, relative_locus)
    unmapped_loci = read_kmer(contig_interval)
    print (unmapped_loci)
    my_breakpoints = get_breakpoints(unmapped_loci)
    print (my_breakpoints)
    split_fasta(input_fasta, output_fasta, my_breakpoints)
