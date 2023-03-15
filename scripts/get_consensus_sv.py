"""
Generate consensus fasta file from results of PBSV2
"""
from Bio import SeqIO
import vcf
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys

# ref = "/mnt/d/breakpoints/assembly/sim/database/ecoli/NZ_CP049085.2.fasta"
# consensus_file = "/home/wangshuai/assembly_result/pacbio/test.fasta"
# sv_vcf = "/home/wangshuai/assembly_result/pacbio/ref.var.vcf"

ref= sys.argv[1]
sv_vcf= sys.argv[2]
consensus_file = sys.argv[3]



def load_origin_ref():
    origin_dict = {}
    with open(ref, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            origin_dict[record.id] = str(record.seq)
            # print(record.id)
            # print(record.seq)
    return origin_dict



def read_vcf():
    new_dict = origin_dict.copy()
    with open(sv_vcf, "r") as vcf_file:
        vcf_reader = vcf.Reader(vcf_file)

        for record in vcf_reader:
            if not record.is_sv:
                continue
            if record.INFO['SVTYPE'] == "BND":
                continue

            if record.INFO['SVTYPE'] == "DEL":
                start = len(origin_dict[record.CHROM]) - record.POS - 1
                end = start + abs(int(record.INFO['SVLEN'][0]))
                start *= -1
                end *= -1
                new_dict[record.CHROM] = new_dict[record.CHROM][:start] + new_dict[record.CHROM][end:]
            
            elif record.INFO['SVTYPE'] == "INS":
                start = len(origin_dict[record.CHROM]) - record.POS - 1
                start *= -1
                new_dict[record.CHROM] = new_dict[record.CHROM][:start] + str(record.ALT[0]) + new_dict[record.CHROM][start:]

            elif record.INFO['SVTYPE'] == "INV":
                start = len(origin_dict[record.CHROM]) - record.POS - 1
                end = start + abs(int(record.INFO['SVLEN'][0]))
                start *= -1
                end *= -1
                reverse_complement = str(Seq(new_dict[record.CHROM][start:end]).reverse_complement())
                new_dict[record.CHROM] = new_dict[record.CHROM][:start] + reverse_complement + new_dict[record.CHROM][end:]

            print(record.CHROM, record.POS, record.INFO['SVTYPE'], record.INFO['SVLEN'], record.ALT)
    
    return new_dict

def output():
    new_records = []

    for chrom in new_dict:
        record = SeqRecord(Seq(new_dict[chrom]), chrom, "", "")
        new_records.append(record)
    SeqIO.write(new_records, consensus_file, "fasta")

origin_dict = load_origin_ref()
new_dict = read_vcf()
output()