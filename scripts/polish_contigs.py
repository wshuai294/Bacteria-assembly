import sys, os, re
from Bio import SeqIO

min_sv_len = 50

def read_delly_vcf(vcffile):
    my_breakpoints = []
    remove_repeat = {}
    if not os.path.isfile(vcffile):
        print ("no sv vcf.")
        return 0
    for line in open(vcffile, "r"):
        line = line.strip()
        if line[0] == "#":
            continue

        if re.search("IMPRECISE", line):
            continue
        end_pos=re.search(";END=(.*?);", line)
        array = line.split()
        chrom = array[0]
        pos = int(array[1])

        if not re.search("SVTYPE=BND", line):
            sv_len = abs(int(end_pos.group(1)) - pos)
            if sv_len < min_sv_len:
                continue
            if chrom + "_" + str(pos) not in remove_repeat:
                my_breakpoints.append([chrom, pos])
                remove_repeat[chrom + "_" + str(pos)] = 1
            if array[4] == "<DEL>" or array[4] == "<DUP>":
                pos_end = pos + sv_len
            if chrom + "_" + str(pos) not in remove_repeat:
                my_breakpoints.append([chrom, pos_end])
                remove_repeat[chrom + "_" + str(pos_end)] = 1
        else:
            fir = re.search("\[(.*?)\[", array[4])
            sec = re.search("](.*?)]", array[4])
            if fir:
                other_end = fir.group(1)
            else:
                other_end = sec.group(1)
            locus_info = other_end.split(":")
            chrom2 = locus_info[0]
            pos2 = int (locus_info[1])
            if chrom2 == chrom and abs(pos2 - pos) < min_sv_len:
                continue
            if chrom + "_" + str(pos) not in remove_repeat:
                my_breakpoints.append([chrom, pos])
                remove_repeat[chrom + "_" + str(pos)] = 1
            if chrom2 + "_" + str(pos2) not in remove_repeat:
                my_breakpoints.append([chrom2, pos2])
                remove_repeat[chrom2 + "_" + str(pos2)] = 1
    return sorted(my_breakpoints)

def split_fasta(input_fasta, output_fasta, breakpoint_list):
    # Load the FASTA file
    fasta_file = SeqIO.parse(input_fasta, "fasta")
    handle = open(output_fasta, "w")
    # Split the FASTA file based on the breakpoints
    for record in fasta_file:
        start = 0
        end = len(record)
        for breakpoint in breakpoint_list:
            if breakpoint[0] == record.id:
                if abs(start - breakpoint[1]) > min_sv_len:
                    new_record = record[start:breakpoint[1]]
                    new_record.id = f"{record.id}_{start}_{breakpoint[1]}"
                    new_record.description = ""
                    SeqIO.write(new_record, handle, "fasta")
                    start = breakpoint[1]
        if abs(start - end) > min_sv_len:
            new_record = record[start:end]
            new_record.id = f"{record.id}_{start}_{end}"
            new_record.description = ""
            SeqIO.write(new_record, handle, "fasta")
    
vcffile = sys.argv[1]
input_fasta = sys.argv[2]
output_fasta = sys.argv[3]
my_breakpoints = read_delly_vcf(vcffile)
# print (my_breakpoints)
split_fasta(input_fasta, output_fasta, my_breakpoints)