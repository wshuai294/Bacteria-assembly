"""
filter assembled contigs according to length
"""

from Bio import SeqIO
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import os


original_assemblies = sys.argv[1]
novel_assemblies = sys.argv[2]
min_length = int(sys.argv[3])


def main():
    
    long_sequences = []  # Setup an empty list
    for record in SeqIO.parse(original_assemblies, "fasta"):
        if len(record.seq) > min_length:
            # Add this record to our list
            long_sequences.append(record)

    print("Found %i short sequences" % len(long_sequences))

    SeqIO.write(long_sequences, novel_assemblies, "fasta")

if os.path.isfile(original_assemblies):
    main()
else:
    print ("No assembled insertions.")