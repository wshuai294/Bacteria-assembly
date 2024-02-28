"""
Given a genome, map reads to it, and assemble the unmapped reads. 

Spades must be installed in the system path

"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import sys
import argparse
import numpy as np
import csv

from module_alignment import alignment


def run():
    bam_prefix = f"{options.o}/{options.s}.split" 
    alignment(options.align_tool, options.genome, options.fq1, options.fq2, bam_prefix, options.t)


    command = f"""
    dir={sys.path[0]}
    fq1={options.fq1}
    fq2={options.fq2}
    split_ref={options.genome}
    threads={options.t}
    sample={options.o}/{options.s}
    outdir={options.o}
    seg_ref={options.o}/{options.s}.segment.fasta
    assembly_dir={options.o}/{options.s}_assembly

    # samtools view -f 4 -b $sample.split.bam > $sample.unmap.bam
    python $dir/extract_unmap.py $sample.split.bam $sample.unmap.bam
    samtools fastq -1 $sample.unmapped.1.fq -2 $sample.unmapped.2.fq -s $sample.unmapped.s.fq $sample.unmap.bam
    gzip -f $sample.unmapped.*fq
    
    cat $split_ref>$seg_ref

    rm -r $assembly_dir
    echo "start spades..."
    spades.py --isolate -t $threads -1 $sample.unmapped.1.fq.gz -2 $sample.unmapped.2.fq.gz -s $sample.unmapped.s.fq.gz -o $assembly_dir >$sample.spades.log

    if [ -f "$assembly_dir/contigs.fasta" ]; then
        # cat $assembly_dir/contigs.fasta >>$seg_ref
        python $dir/filter_assemblies.py $assembly_dir/contigs.fasta $assembly_dir/contigs.filter.fasta {options.m} 10
        cat $assembly_dir/contigs.filter.fasta >>$seg_ref
        # contig=$assembly_dir/contigs.filter.fasta
        # $dir/../tools/novoindex $contig.ndx $contig
        # samtools faidx $contig
        # $dir/../tools/novoalign -d $contig.ndx -f $fq1 $fq2 -F STDFQ -o SAM \
        #     -o FullNW --mCPU $threads | samtools view -Sb - | samtools sort -  > $sample.contig.bam
        # samtools index $sample.contig.bam
        # python $dir/find_INS_breakpoints.py $sample.contig.bam $contig $assembly_dir/contigs.splitted.fasta
        # cat $assembly_dir/contigs.splitted.fasta >>$seg_ref
    else
        echo "no assembled contigs"
    fi

    """
    os.system(command)
    # filter_contig = f"{options.o}/{options.s}_assembly/contigs.filter.fasta"
    # if os.path.isfile(filter_contig):
    #     bam_prefix = f"{options.o}/{options.s}.contig" 
    #     alignment(options.align_tool, filter_contig, options.fq1, options.fq2, bam_prefix, options.t)
    #     cmd = f"""
    #     dir={sys.path[0]}
    #     assembly_dir={options.o}/{options.s}_assembly
    #     contig=$assembly_dir/contigs.filter.fasta
    #     seg_ref={options.o}/{options.s}.segment.fasta
    #     sample={options.o}/{options.s}

    #     python $dir/find_INS_breakpoints.py $sample.contig.bam $contig $assembly_dir/contigs.splitted.fasta
    #     cat $assembly_dir/contigs.splitted.fasta >>$seg_ref
    #     """
    #     os.system(cmd)
    # else:
    #     print ("no assembled contigs")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get unmap or repeat segments.", add_help=False, \
    usage="%(prog)s -h", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")
    required.add_argument("--genome", type=str, help="<str> the fasta file.", metavar="\b")
    required.add_argument("--fq1", type=str, help="<str> unzipped fastq 1 file.", metavar="\b")
    required.add_argument("--fq2", type=str, help="<str> unzipped fastq 2 file.", metavar="\b")
    required.add_argument("-s", type=str, default="sample", help="<str> Sample name.", metavar="\b")
    required.add_argument("-o", type=str, default="./", help="<str> Output folder.", metavar="\b")
    optional.add_argument("-m", type=int, default=1000, help="<int> minimum segment length.", metavar="\b")
    optional.add_argument("-t", type=int, default=10, help="<int> number of threads.", metavar="\b")
    optional.add_argument("--align_tool", type=str, default="novoalign", help="alignment method, novoalign or bwa.", metavar="\b")
    optional.add_argument("-h", "--help", action="help")

    options = parser.parse_args()

    if len(sys.argv) == 1:
        print (f"see python {sys.argv[0]} -h")
    else:
        run()