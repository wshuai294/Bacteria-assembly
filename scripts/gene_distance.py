import re
from collections import defaultdict
from Bio import SeqIO
import sys
import os

class Gene(object):

    def __init__(self):
        self.contig = ''
        self.pos = ''
        self.name = 'n/a'
        self.direc = ''
    
    def get_gene(self, gff_line):
        array = gff_line.strip().split()
        # print (array)
        self.contig = array[0]
        self.pos = round((int(array[3]) + int(array[4]))/2)
        self.direc = array[6]
        mat = re.search(";gene=(.*?);inference", gff_line)
        if mat:
            self.name = mat.group(1)

def read_gff(ref_gff):
    gene_info_dict = {} # gene_name: gene_object
    gene_copy_dict = defaultdict(int)  # gene_name: gene copy number in the genome
    contig_gene_dict = defaultdict(list) # contigs_name: the name of genes in that contig

    f = open(ref_gff, 'r')
    for gff_line in f:

        if gff_line.strip() == "##FASTA":
            break
        if gff_line[0] == "#":
            continue

        gene = Gene()
        gene.get_gene(gff_line)

        gene_info_dict[gene.name] = gene
        gene_copy_dict[gene.name] += 1
        contig_gene_dict[gene.contig].append(gene.name)
        # print (gene.name)
    f.close()

    return gene_info_dict, gene_copy_dict, contig_gene_dict

def get_seg_distance(gene_info_dict, gene_copy_dict, contig_gene_dict, gene_info_dict_seg, gene_copy_dict_seg, contig_gene_dict_seg, conver_name_dict, gene_graph):
    h = open(gene_graph, 'w')
    seg_1 = "gnl|X|HEOBEEFN_1"
    seg_2 = "gnl|X|HEOBEEFN_3"
    seg_list = list(contig_gene_dict_seg.keys())
    for i in range(len(seg_list)):
        for j in range(i+1, len(seg_list)):
            seg_1 = seg_list[i]
            seg_2 = seg_list[j]

            # quantify the distance between two contigs
            min_distance = float('inf')
            for gene_name_1 in contig_gene_dict_seg[seg_1]:
                for gene_name_2 in contig_gene_dict_seg[seg_2]:
                    if gene_name_1 == gene_name_2:
                        continue
                    if gene_name_1 not in gene_info_dict or gene_name_2 not in gene_info_dict:
                        continue
                    if gene_copy_dict[gene_name_1] > 1 or gene_copy_dict[gene_name_2] > 1:
                        continue
                    if gene_copy_dict_seg[gene_name_1] > 1 or gene_copy_dict_seg[gene_name_2] > 1:
                        continue
                    if gene_info_dict[gene_name_1].contig != gene_info_dict[gene_name_2].contig:
                        print (f"{gene_name_1} and {gene_name_2} locate in different contigs on the ref.")
                    ## cal gene distance on ref
                    distance = abs(gene_info_dict[gene_name_1].pos - gene_info_dict[gene_name_2].pos)
                    if distance < min_distance:
                        min_distance = distance
            if min_distance != float('inf'):
                print ("JUNC", conver_name_dict[seg_1], "+", conver_name_dict[seg_2], "+", min_distance, file = h)
    h.close()

def get_contig_names(fasta_file):
    contig_names = []

    # Read the FASTA file and extract contig names
    for record in SeqIO.parse(fasta_file, "fasta"):
        contig_names.append(record.id)

    return contig_names

def rename_contigs(fasta_file, prefix, rename_seg):
    renamed_records = []

    # Read the FASTA file and rename contigs
    i = 1
    for record in SeqIO.parse(fasta_file, "fasta"):
        new_name = f"{prefix}_{i}"  # Modify the prefix as desired
        record.id = new_name
        record.description = ""
        renamed_records.append(record)
        i += 1

    # Save the renamed contigs in a new FASTA file
    output_file = rename_seg
    with open(output_file, "w") as f:
        SeqIO.write(renamed_records, f, "fasta")

    print("Contigs have been renamed and saved to 'renamed.fasta'.")


def rename(ori_seg, rename_seg):
    conver_name_dict = {}
    ori_name = get_contig_names(ori_seg)
    new_name = get_contig_names(rename_seg)
    # print (len(new_name), len(ori_name))
    for i in range(len(ori_name)):
        conver_name_dict[new_name[i]] = ori_name[i]
    return conver_name_dict

def run_prokka(fasta, outdir, prefix):
    cmd = f"""prokka --quiet --force --cpus {cpus} --outdir {outdir}  --prefix {prefix} {fasta}"""
    os.system(cmd)


if __name__ == "__main__":    

    ref = sys.argv[1]
    outdir = sys.argv[2]
    seg_fa = sys.argv[3]
    ID = sys.argv[4]
    cpus = sys.argv[5]
    rename_seg = seg_fa + ".rename.fasta"
    gene_graph = outdir + "/" + ID + ".gene.graph.txt"

    if os.path.isfile(ref + ".gff"):
        ref_gff = ref + ".gff"
    else:
        
        ref_gff = outdir + "/" + ID + "_ref.gff" 
        if  not os.path.isfile(ref_gff):
            run_prokka(ref, outdir, ID+"_ref")
        print (ref_gff)


    # ref_gff = "/home/wangshuai/assembly_result/truth/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2.gff"
    # seg_gff = "/home/wangshuai/assembly_result/w6_ecoli/seg/seg.gff"

    # ori_seg = "/home/wangshuai/assembly_result/w6_ecoli/NZ_AP022525.1.seg.fasta"
    # rename_seg = "/home/wangshuai/assembly_result/w6_ecoli/NZ_AP022525.1.seg.rename.fasta"

    # gene_graph = "/home/wangshuai/assembly_result/w6_ecoli/gene.graph.txt"

    rename_contigs(seg_fa, "seg", rename_seg)
    conver_name_dict = rename(seg_fa, rename_seg)
    run_prokka(rename_seg, outdir, ID+"_seg")
    seg_gff = outdir + "/" + ID + "_seg.gff" 

    gene_info_dict, gene_copy_dict, contig_gene_dict = read_gff(ref_gff)
    gene_info_dict_seg, gene_copy_dict_seg, contig_gene_dict_seg = read_gff(seg_gff)
    # print (contig_gene_dict_seg.keys())
    get_seg_distance(gene_info_dict, gene_copy_dict, contig_gene_dict, gene_info_dict_seg, gene_copy_dict_seg, contig_gene_dict_seg, conver_name_dict, gene_graph)


    