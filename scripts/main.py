import os
import sys
import argparse

def check_input():

    if not os.path.exists(options.o):
        os.system(f"mkdir {options.o}")
    if not os.path.exists(options.genome_list):
        print ("Reference list file: %s is not detected."%(options.genome_list))
        sys.exit()
    if not os.path.exists(options.fq1):
        print ("Reference file: %s is not detected."%(options.fq1))
        sys.exit()


parser = argparse.ArgumentParser(description="Assembly.", add_help=False, \
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
# optional.add_argument("--sample_ratio", type=int, default=100, help="<0-100> sample ratio (%).", metavar="\b")
optional.add_argument("-m", type=int, default=1000, help="<int> minimum segment length.", metavar="\b")
optional.add_argument("-p", type=int, default=1, help="<0/1> whether refine the contigs with Pilon.", metavar="\b")
optional.add_argument("--sample_depth", type=float, default=5, help="<float> only retain reads with this depth in downsampling.", metavar="\b")
optional.add_argument("--distance", type=int, default=10, help="<int> merge two interval with distance shorted than this value.", metavar="\b")
optional.add_argument("--align_tool", type=str, default="novoalign", help="alignment method, novoalign or bwa.", metavar="\b")

optional.add_argument("-h", "--help", action="help")

options = parser.parse_args()

if len(sys.argv) == 1:
    print (f"see python {sys.argv[0]} -h")
else:
    check_input()

    cmd = f"""
        python {sys.path[0]}/01_find_best_ref.py --genome_list {options.genome_list} --fq1 {options.fq1} --fq2 {options.fq2} \
        -s {options.s} -o {options.o} -k {options.k} -t {options.t} -d {options.d} --sample_depth {options.sample_depth}
        """
    os.system(cmd)

    ref_cmd = f"sample={options.o}/{options.s}"
    ref_cmd +=  """
    highest=$(awk -F',' 'BEGIN { max = 0 } 
            { if($2>max) { max=$2; val=$1; second_col=$2 } } 
            END { print val "," second_col }' $sample.match_rate.csv)
    ref=$(echo $highest | awk -F',' '{print $1}')
    echo "selected ref is $highest."
    echo "$ref" >$sample.selected.ref.txt
    """
    os.system(ref_cmd)

    cmd = f"""
        best_ref=$(cat {options.o}/{options.s}.selected.ref.txt)

        python {sys.path[0]}/02_find_ref_segments.py -r $best_ref --fq1 {options.fq1} --fq2 {options.fq2} \
        -s {options.s} -o {options.o} -t {options.t} -m {options.m} --distance {options.distance} -k {options.k}

        python {sys.path[0]}/03_get_unmap_segments.py --genome {options.o}/{options.s}.split.fasta --fq1 {options.fq1} --fq2 {options.fq2} \
        -s {options.s} -o {options.o} -t {options.t} --align_tool {options.align_tool} -m {options.m}

        python {sys.path[0]}/04_ref_scaffold.py --ref_genome $best_ref --fq1 {options.fq1} --fq2 {options.fq2} \
        -s {options.s} -o {options.o} -t {options.t} --align_tool {options.align_tool} --seg_ref {options.o}/{options.s}.segment.fasta -p {options.p}

    """
    os.system(cmd)