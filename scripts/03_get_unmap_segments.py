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


def run():
    command = f"""
    dir={sys.path[0]}
    fq1={options.fq1}
    fq2={options.fq2}
    split_ref={options.genome}
    threads={options.t}
    sample={options.o}/{options.s}
    outdir={options.o}
    seg_ref={options.o}/{options.s}.segment.fasta

    $dir/../tools/novoindex $split_ref.ndx $split_ref
    samtools faidx $split_ref
    $dir/../tools/novoalign -d $split_ref.ndx -f $fq1 $fq2 -F STDFQ -o SAM \
        -o FullNW --mCPU $threads | samtools view -Sb - | samtools sort -  > $sample.split.bam
    samtools index $sample.split.bam
    samtools view -f 4 -b $sample.split.bam > $sample.unmap.bam
    samtools fastq -1 $sample.unmapped.1.fq -2 $sample.unmapped.2.fq -s $sample.unmapped.s.fq $sample.unmap.bam
    gzip -f $sample.unmapped.*fq
    rm -r $outdir/ass
    cat $split_ref>$seg_ref


    echo "start spades..."
    spades.py --isolate -t $threads -1 $sample.unmapped.1.fq.gz -2 $sample.unmapped.2.fq.gz -s $sample.unmapped.s.fq.gz -o $outdir/ass >$sample.spades.log

    if [ -f "$outdir/ass/contigs.fasta" ]; then
        python $dir/filter_assemblies.py $outdir/ass/contigs.fasta $outdir/ass/contigs.filter.fasta 100 10
        contig=$outdir/ass/contigs.filter.fasta
        $dir/../tools/novoindex $contig.ndx $contig
        samtools faidx $contig
        $dir/../tools/novoalign -d $contig.ndx -f $fq1 $fq2 -F STDFQ -o SAM \
            -o FullNW --mCPU $threads | samtools view -Sb - | samtools sort -  > $sample.contig.bam
        samtools index $sample.contig.bam
        python $dir/find_INS_breakpoints.py $sample.contig.bam $contig $outdir/ass/contigs.splitted.fasta
        cat $outdir/ass/contigs.splitted.fasta >>$seg_ref
    else
        echo "no assembled contigs"
    fi

    """
    os.system(command)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get unmap or repeat segments.", add_help=False, \
    usage="%(prog)s -h", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")
    required.add_argument("--genome", type=str, help="<str> ba fasta file.", metavar="\b")
    required.add_argument("--fq1", type=str, help="<str> unzipped fastq 1 file.", metavar="\b")
    required.add_argument("--fq2", type=str, help="<str> unzipped fastq 2 file.", metavar="\b")
    required.add_argument("-s", type=str, default="sample", help="<str> Sample name.", metavar="\b")
    required.add_argument("-o", type=str, default="./", help="<str> Output folder.", metavar="\b")

    optional.add_argument("-t", type=int, default=10, help="<int> number of threads.", metavar="\b")
    optional.add_argument("-h", "--help", action="help")

    options = parser.parse_args()

    if len(sys.argv) == 1:
        print (f"see python {sys.argv[0]} -h")
    else:
        run()