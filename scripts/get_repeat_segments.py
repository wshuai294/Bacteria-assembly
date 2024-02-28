import subprocess
from Bio import SeqIO

def self_alignment_blast(genome_file, blast_output, repeat_output):
    # Run BLAST command
    blast_command = f"blastn -query {genome_file} -subject {genome_file} -outfmt 6 -out {blast_output}"
    subprocess.run(blast_command, shell=True)
    
    # Extract repeat segments from BLAST results
    repeat_segments = {}
    
    genome_records = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
    
    with open(blast_output) as result_handle:
        result_handle.readline()
        for line in result_handle:
            field =  line.strip().split("\t")
            start, end = int(field[6]), int(field[7])
            chrom = field[0]
            repeat_segment = genome_records[chrom].seq[start-1:end].upper()
            # Consider repeat segments with length over 1000 bp
            segname = f"repeat_{chrom}_{start}_{end}"
            if len(repeat_segment) > 1000:
                repeat_segments[segname] = repeat_segment
    
    # Save repeat segments in a new FASTA file
    with open(repeat_output, "w") as output_handle:
        for segname in repeat_segments:
            output_handle.write(f">Repeat_{segname}\n{repeat_segments[segname]}\n")

# Example usage
genome_file = "/mnt/d/breakpoints/assembly/sim/database/Escherichia_coli/NZ_CP008957.1.fasta"
blast_output = "/home/wangshuai/assembly_result_v2/test/blast_results.txt"
repeat_output = "/home/wangshuai/assembly_result_v2/test/repeat_segments.fasta"

self_alignment_blast(genome_file, blast_output, repeat_output)
