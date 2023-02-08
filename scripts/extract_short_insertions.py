import pysam
import sys

short_indel_vcf = sys.argv[1]
short_indel_fasta = sys.argv[2]
short_indel_pos = sys.argv[3]

vcf = pysam.VariantFile(short_indel_vcf)
fasta = open(short_indel_fasta, "w")
meta = open(short_indel_pos, "w")

for record in vcf.fetch():
    ref = record.ref
    alt = record.alts[0]

    if len(ref) < len(alt):
        ins_seq = alt[len(ref):]
        if len(ins_seq) > 50:
            fasta.write(">NODE_INS_" + str(record.pos) + "_" + str(len(ins_seq)) +"_cov_xx" + "\n" + ins_seq + "\n")
            meta.write(record.chrom + "\t" + str(record.pos) + "\n")

vcf.close()
fasta.close()
meta.close()
