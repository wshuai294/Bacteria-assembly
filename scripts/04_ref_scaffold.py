import sys
import math
import pysam
import numpy as np
import random
import re
import os
import argparse
import time
from Bio import SeqIO
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib

from module_alignment import alignment
from module_gene_graph import get_gene_graph

def run_ragtag():
    cmd = f"""
    outdir={options.o}
    sample={options.s}
    ref={options.ref_genome}
    # ragtag.py correct -t {options.t} -u -o $outdir/ragtag_$sample $ref $outdir/$sample.segment.fasta
    # ragtag.py scaffold -t {options.t} -w -u -o $outdir/ragtag_$sample $ref $outdir/ragtag_$sample/ragtag.correct.fasta

    ragtag.py scaffold -t {options.t} -w -u -o $outdir/ragtag_$sample $ref $outdir/$sample.segment.fasta
    """ 
    os.system(cmd)

def refine_assembly():
    raw_contig = f"{options.o}/{options.s}.contigs.fasta"
    bam_prefix = f"{options.o}/{options.s}.contigs"
    alignment(options.align_tool, raw_contig, options.fq1, options.fq2, bam_prefix, options.t)
    cmd = f"""
    sample={options.o}/{options.s}
    samtools index $sample.contigs.bam
    pilon -Xmx10g --genome $sample.contigs.fasta --frags $sample.contigs.bam --output $sample.contigs.refine --chunksize 1000000
    """
    os.system(cmd)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Scaffold with reference.", add_help=False, \
    usage="%(prog)s -h", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")
    required.add_argument("--seg_ref", type=str, help="<str> segment reference, a fasta file.", metavar="\b")
    required.add_argument("--fq1", type=str, help="<str> unzipped fastq 1 file.", metavar="\b")
    required.add_argument("--fq2", type=str, help="<str> unzipped fastq 2 file.", metavar="\b")
    required.add_argument("-s", type=str, default="sample", help="<str> Sample name.", metavar="\b")
    required.add_argument("-o", type=str, default="./", help="<str> Output folder.", metavar="\b")
    required.add_argument("--ref_genome", type=str, help="<str> reference genome to get gene graph, fasta file.", metavar="\b")

    optional.add_argument("-t", type=int, default=10, help="<int> number of threads.", metavar="\b")
    optional.add_argument("-p", type=int, default=1, help="<0/1> whether refine the contigs with Pilon.", metavar="\b")
    optional.add_argument("--align_tool", type=str, default="novoalign", help="alignment method, novoalign or bwa.", metavar="\b")
    optional.add_argument("-h", "--help", action="help")

    options = parser.parse_args()

    if len(sys.argv) == 1:
        print (f"see python {sys.argv[0]} -h")
    else:

        run_ragtag()
        rag_tag_fasta = f"{options.o}/ragtag_{options.s}/ragtag.scaffold.fasta"
        raw_contig = options.o + "/" + options.s + ".contigs.fasta"
        os.system(f"cp {rag_tag_fasta} {raw_contig}")

        print ("assembly result is in %s"%(raw_contig))
        if options.p == 1:
            print ("refine the contigs...")
            refine_assembly()
            print ("refined assembly result is in %s"%(options.o + "/" + options.s + ".contigs.refine.fasta"))

            