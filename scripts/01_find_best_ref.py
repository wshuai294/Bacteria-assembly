"""
Given a genome list, find the best-match genome as the reference for downstream assembly. 

"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import sys
import argparse
import numpy as np
import csv

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

def estimate_sample_ratio(fasta_file, fastq_file, depth):
    # best depth 16x
    est_depth = float(count_bases_fastq(fastq_file))/count_bases_fasta(fasta_file)
    best_sample_ratio = round(depth/est_depth * 100/2)
    return best_sample_ratio

def run():
    command = f"""
    {sys.path[0]}/src/select_ref {options.fq1} {options.fq2} {options.genome_list} {options.k} {options.d} {options.t} {options.o}/{options.s}.match_rate.csv {best_sample_ratio}  


    """
    os.system(command)
    get_line_with_largest_value()

def get_line_with_largest_value():
    max_value = None
    max_line = None
    csv_file = f"{options.o}/{options.s}.match_rate.csv"
    with open(csv_file, 'r') as file:
        csv_reader = csv.reader(file)
        next(csv_reader)  # Skip header row

        for line in file:
            values = line.strip().split(',')
            value = float(values[1])  # Assuming the second column contains numeric values

            if max_value is None or value > max_value:
                max_value = value
                max_line = line.strip()
                
    best_file = open(f"{options.o}/{options.s}.selected.ref.txt", 'w')
    print (max_line, end = '', file = best_file)
    best_file.close()
    # return max_line


if __name__ == "__main__":


    parser = argparse.ArgumentParser(description="Find the best reference genome.", add_help=False, \
    usage="%(prog)s -h", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")
    required.add_argument("--genome_list", type=str, help="<str> ba fasta file list, where each fasta is a single chromosome", metavar="\b")
    required.add_argument("--fq1", type=str, help="<str> unzipped fastq 1 file.", metavar="\b")
    required.add_argument("--fq2", type=str, help="<str> unzipped fastq 2 file.", metavar="\b")
    required.add_argument("-s", type=str, default="sample", help="<str> Sample name.", metavar="\b")
    required.add_argument("-o", type=str, default="./", help="<str> Output folder.", metavar="\b")
    optional.add_argument("-k", type=int, default=26, help="<int> kmer length.", metavar="\b")
    optional.add_argument("-t", type=int, default=10, help="<int> number of threads.", metavar="\b")
    optional.add_argument("-d", type=int, default=3, help="<int> lowest kmer count that regarded as hit.", metavar="\b")
    optional.add_argument("--sample_depth", type=float, default=5, help="<float> only retain reads with this depth in downsampling.", metavar="\b")
    # optional.add_argument("--sample_ratio", type=int, default=100, help="<0-100> sample ratio (%).", metavar="\b")
    optional.add_argument("-h", "--help", action="help")

    options = parser.parse_args()

    if len(sys.argv) == 1:
        print (f"see python {sys.argv[0]} -h")
    else:
        infile = open(options.genome_list, 'r')
        first_genome = infile.readline().strip()
        best_sample_ratio = estimate_sample_ratio(first_genome, options.fq1, options.sample_depth)
        print (f"sampling ratio is {best_sample_ratio}.")
        run()

        # 