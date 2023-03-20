"""
require jellyfish
"""

import subprocess

# Define the input FASTQ file and the k-mer length
# fastq_file = "~/assembly_result/test_short.1.fq"
fastq_file = ""
kmer_length = "13"

# Count k-mers using Jellyfish
jellyfish_cmd = ["jellyfish", "count", "-m", kmer_length, "-s", "100M", "-t", "8", "-C", fastq_file]
subprocess.run(jellyfish_cmd)

# Estimate sequencing depth using Jellyfish
jellyfish_stats_cmd = ["jellyfish", "stats", "mer_counts.jf"]
jellyfish_process = subprocess.Popen(jellyfish_stats_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
jellyfish_output, jellyfish_error = jellyfish_process.communicate()

# Extract the total number of k-mers and unique k-mers from the Jellyfish output
total_kmers = int(jellyfish_output.split()[1])
unique_kmers = int(jellyfish_output.split()[3])

# Calculate the sequencing depth
seq_depth = total_kmers / unique_kmers

# Print the estimated sequencing depth
print("Estimated sequencing depth:", seq_depth)
