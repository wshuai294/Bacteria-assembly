import sys, os, re
from Bio import SeqIO
import pysam
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

min_sv_len = 50
min_variant_q = 20

def read_delly_vcf(vcffile):
    my_breakpoints = {}
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

        if chrom not in my_breakpoints:
            my_breakpoints[chrom] = []

        if not re.search("SVTYPE=BND", line):
            sv_len = abs(int(end_pos.group(1)) - pos)
            if sv_len < min_sv_len:
                continue
            my_breakpoints[chrom].append(pos)

            if array[4] == "<DEL>" or array[4] == "<DUP>":
                pos_end = pos + sv_len

            my_breakpoints[chrom].append(pos_end)

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

            my_breakpoints[chrom].append(pos)
            if chrom2 not in my_breakpoints:
                my_breakpoints[chrom2] = []
            my_breakpoints[chrom2].append(pos2)
    
    for chrom in my_breakpoints:
        my_breakpoints[chrom] = list(set(my_breakpoints[chrom]))
        my_breakpoints[chrom].sort()
        # print (chrom, my_breakpoints[chrom])
    return my_breakpoints

def split_fasta(input_fasta, output_fasta, my_breakpoints):
    # Load the FASTA file
    fasta_file = SeqIO.parse(input_fasta, "fasta")
    handle = open(output_fasta, "w")
    # Split the FASTA file based on the breakpoints
    for record in fasta_file:
        start = 0
        end = len(record)
        if record.id in my_breakpoints:
            breakpoint_list = my_breakpoints[record.id]
        else:
            breakpoint_list = []
        for breakpoint in breakpoint_list:
            if abs(start - breakpoint) > min_sv_len:
                new_record = record[start:breakpoint]
                new_record.id = f"{record.id}_{start}_{breakpoint}"
                new_record.description = ""
                SeqIO.write(new_record, handle, "fasta")
                start = breakpoint
        if abs(start - end) > min_sv_len:
            new_record = record[start:end]
            new_record.id = f"{record.id}_{start}_{end}"
            new_record.description = ""
            SeqIO.write(new_record, handle, "fasta")

def read_mismatch(vcf, method):

    vcf = pysam.VariantFile(vcf)

    for record in vcf.fetch():
        if record.qual < min_variant_q:
            continue
        ref = record.ref
        alt = record.alts[0]

        if method == "svaba" and abs(len(ref) - len(alt)) < min_sv_len:
            continue
        if method == "freebayes" and abs(len(ref) - len(alt)) >= min_sv_len:
            continue
        if record.chrom not in my_mismatches:
            my_mismatches[record.chrom] = []
        my_mismatches[record.chrom].append([record.chrom, record.pos, ref, alt])
    
    vcf.close()

def get_consensus():
    fasta_file = SeqIO.parse(input_fasta, "fasta")
    handle = open(middle_fasta, "w")

    position_correction = {}

    for record in fasta_file:
        # Create a list to store the consensus sequence
        consensus = []
        position_correction[record.id] = {}
        relative_pos = 0 

        if record.id in my_mismatches:
            variant_list = my_mismatches[record.id]
            variant_list.sort(key = lambda x: x[1])
            original_length = len(record.seq)
            reference_seq = str(record.seq)
            
            # Initialize a start position and end position for each variant
            start = 0
            end = 0
            for chrom, pos, ref, alt in variant_list:
                # Update the end position to be the start of the current variant
                end = pos - 1
                
                # Add the reference sequence between the start and end positions to the consensus
                consensus.append(reference_seq[start:end])
                
                # Update the start position to be the end of the current variant
                start = end + len(ref)
                
                # Add the alternate sequence to the consensus
                consensus.append(alt)

                relative_pos = relative_pos + (len(alt) - len(ref))
                position_correction[record.id][pos] = relative_pos
            
            # Add the remaining reference sequence to the consensus
            consensus.append(reference_seq[start:])
            
            # Join the consensus sequence into a single string
            consensus_seq = "".join(consensus)
            record.seq = Seq(consensus_seq)
        
        SeqIO.write(record, handle, "fasta")
    
    handle.close()
    return position_correction

def correct_sv_position(my_breakpoints, position_correction):
    revised_breakpoints = {}
    for chrom in my_breakpoints:
        revised_breakpoints[chrom] = []

        if chrom not in position_correction:
            revised_breakpoints[chrom] = my_breakpoints[chrom]
            continue

        for original_pos in my_breakpoints[chrom]:
            pos_diff = 0
            for variant_locus in position_correction[chrom]:
                if variant_locus > original_pos:
                    break
                if original_pos >= variant_locus:
                    pos_diff = position_correction[chrom][variant_locus]
            revised_pos = original_pos + pos_diff

            revised_breakpoints[chrom].append(revised_pos)
    return revised_breakpoints



if __name__ == "__main__":

    freebayes_vcf = sys.argv[1]
    svaba_indel = sys.argv[2]
    svaba_sv = sys.argv[3]

    input_fasta = sys.argv[4]
    middle_fasta = input_fasta + ".tmp.fasta"
    output_fasta = sys.argv[5]


    my_mismatches = {}
    read_mismatch(freebayes_vcf, "freebayes")
    read_mismatch(svaba_indel, "svaba")
    my_breakpoints = read_delly_vcf(svaba_sv)
    position_correction = get_consensus()
    # print (my_breakpoints)
    my_breakpoints = correct_sv_position(my_breakpoints, position_correction)
    # print (my_breakpoints)


    split_fasta(middle_fasta, output_fasta, my_breakpoints)