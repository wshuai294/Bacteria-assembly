from Bio import SeqIO
from Bio.Blast import NCBIXML
import subprocess

def remove_redundant_overlap(fasta_file, output_file):
    # Read the FASTA file and store contig sequences
    contigs = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        contigs[record.id] = record.seq

    # Perform BLAST alignment
    blast_output = "/home/wangshuai/assembly_result/w6_ecoli/output.xml"
    blast_cmd = f"blastn -subject {fasta_file} -query {fasta_file} -out {blast_output} -outfmt 5"
    subprocess.run(blast_cmd, shell=True, check=True)

    # Parse BLAST output XML
    with open(blast_output, "r") as blast_file:
        blast_records = NCBIXML.parse(blast_file)
        for blast_record in blast_records:
            query_id = blast_record.query
            for alignment in blast_record.alignments:
                hit_id = alignment.hit_id
                if hit_id != query_id:
                    hit_start = alignment.hsps[0].sbjct_start
                    hit_end = alignment.hsps[0].sbjct_end
                    contigs[hit_id] = contigs[hit_id][hit_start-1:]

    # Write updated contigs to a new FASTA file
    with open(output_file, "w") as file:
        for contig_id, sequence in contigs.items():
            file.write(f">{contig_id}\n")
            file.write(f"{sequence}\n")

# Example usage
fasta_file = "/home/wangshuai/assembly_result/w6_ecoli/NZ_AP022525.1.seg.fasta"
output_file = "/home/wangshuai/assembly_result/w6_ecoli/output.fasta"
remove_redundant_overlap(fasta_file, output_file)