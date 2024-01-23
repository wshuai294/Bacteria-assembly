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

def readFilter(read):
    return (read.is_proper_pair and
            read.is_paired and
            read.tlen > 0 and
            read.tlen < 1000 and
            not read.is_supplementary and
            not read.is_duplicate and
            not read.is_unmapped and
            not read.mate_is_unmapped)

def getInsertSize(unique_bamfile):
    read_length_list = []
    insert_size_list = []
    r_num = 0
    for read in unique_bamfile:
        if readFilter(read):
            insert_size_list.append(read.tlen)
            read_length_list.append(len(read.query_sequence))
            r_num += 1
        if r_num > 10000:
            print ("consider 10000 reads in estimating read length and insert size.")
            break
        
    read_length = int(sum(read_length_list) / len(read_length_list))
    mean = float(sum(insert_size_list)) / len(insert_size_list)
    sdev = math.sqrt(float(sum([(x - mean)**2 for x in insert_size_list])) / (len(insert_size_list) - 1))
    return mean, sdev, read_length, r_num

def map_ratio(read):
    # cal the mapped ratio in a read
    mapped_len = 0
    read_len = 0
    for ci in read.cigar: # record left clipped length
        if ci[0] == 0:
            mapped_len += int(ci[1])
        read_len += int(ci[1])
    return float(mapped_len)/read_len

def rename_contigs(fasta_file):
    contig_len_dict = {}

    # Read the FASTA file and rename contigs
    for record in SeqIO.parse(fasta_file, "fasta"):
        contig_len_dict[record.id] = len(str(record.seq))
    return contig_len_dict

def get_ref_len(ref_name):
    # if len(ref_name.split(':')) < 2:
    #     ref_len = int(ref_name.split('_')[3])
    # elif re.search("%", ref_name):
    #     # print (ref_name, ref_name.split('%'))
    #     ref_len = int(ref_name.split('%')[2]) -  int(ref_name.split('%')[1])
    # else:
    #     ref_len = int(ref_name.split(':')[1].split('-')[1]) - int(ref_name.split(':')[1].split('-')[0])
    # return ref_len
    return contig_len_dict[ref_name]

def calCrossReads(bam_name, graph, read_type="split"): # read_type: split or pair
    # good_read_dict = filter_bam(bam_name)
    f = open(graph, "w")
    h = open(graph+".read.txt", "w")
    edge_dict = {}
    read_list = []
    bamfile = pysam.AlignmentFile(filename = bam_name, mode = 'rb')
    mean, sdev, rlen, r_num = getInsertSize(bamfile)
    rlen = int(rlen)
    insert_size = int(mean + 2*sdev) + 2 * rlen

    maximum_len = rlen  # only onsider reads mapped on the ends
    for read in bamfile.fetch():
        if read.is_unmapped  or read.mate_is_unmapped:
            continue 
        if read.mapping_quality < min_q:
            continue
        if read.has_tag('XA'):
            continue
        # if map_ratio(read) < min_map_ratio:
        #     continue

        if read_type == "split": # only retain split reads, with SA tag
            if not read.has_tag('SA'):
                continue
            else:
                SA_tag = read.get_tag('SA')
                SA_array = SA_tag.split(";")
                if len(SA_array) > 2:
                    continue
                sa_map = SA_array[0]
                array = sa_map.split(",")
                # print (SA_tag, sa_map)
                other_ref = array[0]
                other_pos = int(array[1])
                other_strand = array[2]
                mate_ref = other_ref + " " + other_strand

        elif read_type == "pair":
            if read.has_tag('SA'):  # not considering split reads
                continue
            other_ref = read.next_reference_name
            if not read.mate_is_reverse:
                mate_ref = read.next_reference_name + " -"
            else:
                mate_ref = read.next_reference_name + " +"
        else:
            print ("wrong read type, shoud be split or pair")
            sys.exit()

            # print (sa_map, other_ref, other_pos)

        if (read.reference_name == read.next_reference_name):
            continue

        ref_len = get_ref_len(read.reference_name)
        mate_len = get_ref_len(other_ref)

        if read_type == "split":  # the read should map to the segment ends, distance with end shorter than read length

            if abs(read.reference_start) > maximum_len and abs(ref_len - read.reference_start) > maximum_len:
                # print (read.query_name, read.reference_start, abs(ref_len - read.reference_start), read.reference_start, ref_len)
                continue
            if not (abs(other_pos) < maximum_len or abs(mate_len - other_pos)< maximum_len):
                continue

        elif read_type == "pair":  # the read should map to the segment ends, distance with end shorter than insert size
            if not (abs(read.reference_start) < insert_size or abs(ref_len - read.reference_start)< insert_size):
                continue
            if not (abs(read.next_reference_start) < insert_size or abs(mate_len - read.next_reference_start)< insert_size):
                continue


        if read.reference_name not in chrom_copy or other_ref not in chrom_copy:
            print ("WARNING: JUNC not in SEG!!!", read.reference_name, other_ref)
            sys.exit()
        if  chrom_copy[read.reference_name] == 0 or  chrom_copy[other_ref] == 0:
            # print ("The copy of one seg is zero in the edge, thus deleted.", read.reference_name, other_ref)
            continue
        # if re.search("NODE_", read.reference_name) and re.search("NODE_", other_ref): # delete edge between two contigs
        #     continue

        if read.is_reverse:
            refname = read.reference_name + " -"
        else:
            refname = read.reference_name + " +"

        edge = refname + " " + mate_ref
        if edge not in edge_dict:
            edge_dict[edge] = 1
        else:
            edge_dict[edge] += 1

        read_list.append([read.query_name, edge])
        
    for chrom in chrom_copy:
        # print ("#\t", chrom, chrom_copy[chrom])
        if chrom_copy[chrom] > 0:
            print ("SEG ", chrom, round(chrom_depth[chrom]), chrom_copy[chrom], 0, 1, file = f)
        else:
            print ("Copy number of %s is zero, thus deleted."%(chrom), round(chrom_depth[chrom]))
    
    # remove_edges, my_nodes = remove_small_circle(edge_dict, chrom_copy, insert_size)
    for edge in edge_dict:
        # if edge_dict[edge] < min_edge_dp:
        #     print ("removed JUNC due to low depth ", edge, round(edge_dict[edge]))
        #     continue
        # if edge in remove_edges:
        #     print ("removed JUNC due to triangle ", edge, round(edge_dict[edge]))
        #     continue
        print ("JUNC ", edge, round(edge_dict[edge]), file = f)
        # print (edge, edge_dict[edge])
    for read_record in read_list:
        print ("READ ", read_record[0], read_record[1], file = h)

    f.close()
    h.close()

def cal_copy_number():
    f = open(depth_file, 'r')
    depth_dict = {}
    all_depth = []
    for line in f:
        array = line.strip().split()
        chrom = array[0]
        pos = int(array[1])
        depth = int(array[2])
        if chrom not in depth_dict:
            depth_dict[chrom] = []
        depth_dict[chrom].append(depth)
        all_depth.append(depth)
    f.close()
    median_depth = np.median(all_depth)
    chrom_copy = {} 
    chrom_depth = {}
    for chrom in depth_dict:
        chrom_depth[chrom] = np.median(depth_dict[chrom])
        chrom_copy[chrom] = round(float(np.median(depth_dict[chrom]))/median_depth)
    return chrom_copy, median_depth, chrom_depth

def remove_small_circle(edge_dict, chrom_copy, insert_size):
    #triangle led by the short middle seg
    seg_length_dict = get_seg_len(chrom_copy)
    # discard
    my_graph = {}
    my_nodes = {}
    for edge in edge_dict:
        if edge_dict[edge] < min_edge_dp:
            continue
        array = edge.split()
        node1 = array[0] + " " + array[1]
        node2 = array[2] + " " + array[3]
        my_nodes[array[0]] = 1
        my_nodes[array[2]] = 1
        if node1 not in my_graph:
            my_graph[node1] = [node2]
        else:
            my_graph[node1] += [node2]
    remove_edges = {}
    for node1 in my_graph:
        if len(my_graph[node1]) == 1:
            pass
        else:
            remove_edge = []
            for node2 in my_graph[node1]:
                if node2 not in my_graph:
                    continue
                else:
                    for node3 in my_graph[node2]:
                        if node3 in my_graph[node1] and seg_length_dict[node2[:-1].strip()] < insert_size:
                            edge = node1 + " " + node3
                            remove_edges[edge] = 1
    return remove_edges, my_nodes

def get_seg_len(chrom_copy):
    # calculate segment length according to its name
    # print (chrom_copy)
    seg_length_dict = {}
    for seg in chrom_copy:
        # array = seg.split("_")
        # if array[0] != "NODE":
        #     new_array = seg[:-1].split(":")[1].split("-")
        #     # print (new_array)
        #     length = abs(int(new_array[0]) - int(new_array[1]))  
        # else:
        #     length = int(array[3])
        # seg_length_dict[seg] = length
        seg_length_dict[seg] = get_ref_len(seg)
    # print (seg_length_dict)
    return seg_length_dict

def plot_graph(graph_file, figure_file):
    g = open(graph_file)
    nodes_list = []
    edges_list = []
    first_partition_nodes = []
    first_partition_segs = {}
    for line in g:  
        array = line.strip().split()
        if array[0] == "SEG":
            continue
        seg1 = array[1]
        node1 = array[1] + array[2]
        seg2 = array[3]
        node2 = array[3] + array[4]
        if seg1 not in first_partition_segs and seg2 not in first_partition_segs:
            first_partition_segs[seg1] = 1
            first_partition_nodes.append(node1)
        elif seg1 in first_partition_segs:
            first_partition_nodes.append(node1)
        elif seg2 in first_partition_segs:
            first_partition_nodes.append(node2)           
        else:
            print ("both in the same side")
        nodes_list += [node1, node2]
        edges_list.append(tuple([node1, node2]))
    first_partition_nodes = sorted(first_partition_nodes)
    DG = nx.Graph()
    # DG = nx.DiGraph()
    DG.add_nodes_from(nodes_list)
    DG.nodes()
    print (len(edges_list), len(list(DG)), len(nodes_list))

    DG.add_edges_from(edges_list, weight = 0.1)
    
    pos = nx.spring_layout(DG, seed=8)    
    # pos = nx.planar_layout(DG) 
    # pos = nx.kamada_kawai_layout(DG) 
    # nx.draw(DG, pos, with_labels=True, node_size=50, width=1, font_size=5)
    options = {
    'node_color': 'blue',
    'node_size': 2,
    'width': 0.2,
    'font_size': 2,
    'arrowstyle': '-|>',
    'arrowsize': 2,
    'alpha': 0.5,
    }
    nx.draw_networkx(DG, pos, arrows=True, **options)

    plt.savefig(figure_file)
    DG.clear()
    g.close()


def run_match(gene_graph_flag):
    ## https://github.com/panguangze/seqGraph/tree/single
    if gene_graph_flag:
        cmd = f"""
            sample={options.o}/{options.s}
            dir={sys.path[0]}

            {options.match_tool} -b --model 1 -v 1 -g $sample.split.graph -r $sample.solve.path.txt \
            -c $sample.solve.c.path.txt -m $sample.new.graph.txt --break_c --gene_junc $sample.gene.graph.txt
            python3 $dir/graph2contig.py {options.seg_ref} $sample.solve.path.txt $sample.contigs.fasta
        """
    else:
        cmd = f"""
            sample={options.o}/{options.s}
            dir={sys.path[0]}

            {options.match_tool} -b --model 1 -v 1 -g $sample.split.graph -r $sample.solve.path.txt \
            -c $sample.solve.c.path.txt -m $sample.new.graph.txt --break_c 
            python3 $dir/graph2contig.py {options.seg_ref} $sample.solve.path.txt $sample.contigs.fasta
        """
    print (cmd)
    os.system(cmd)

def refine_assembly():
    raw_contig = f"{options.o}/{options.s}.contigs.fasta"
    bam_prefix = f"{options.o}/{options.s}.contigs"
    alignment(options.align_tool, raw_contig, options.fq1, options.fq2, bam_prefix, options.t)
    cmd = f"""
    sample={options.o}/{options.s}
    samtools index $sample.contigs.bam
    pilon --genome $sample.contigs.fasta --frags $sample.contigs.bam --output $sample.contigs.refine --chunksize 1000000
    """
    os.system(cmd)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Construct the linkage graph of the segments.", add_help=False, \
    usage="%(prog)s -h", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")
    required.add_argument("--seg_ref", type=str, help="<str> segment reference, a fasta file.", metavar="\b")
    required.add_argument("--fq1", type=str, help="<str> unzipped fastq 1 file.", metavar="\b")
    required.add_argument("--fq2", type=str, help="<str> unzipped fastq 2 file.", metavar="\b")
    required.add_argument("-s", type=str, default="sample", help="<str> Sample name.", metavar="\b")
    required.add_argument("-o", type=str, default="./", help="<str> Output folder.", metavar="\b")

    optional.add_argument("-t", type=int, default=10, help="<int> number of threads.", metavar="\b")
    optional.add_argument("-p", type=int, default=0, help="<0/1> whether refine the contigs with Pilon.", metavar="\b")
    optional.add_argument("-g", type=int, default=0, help="<0/1> whether use the info of the gene distance in reference genome.", metavar="\b")
    required.add_argument("--ref_genome", type=str, help="<str> reference genome to get gene graph, fasta file.", metavar="\b")
    optional.add_argument("--align_tool", type=str, default="novoalign", help="alignment method, novoalign or bwa.", metavar="\b")
    optional.add_argument("--match_tool", type=str, default="/home/wangshuai/softwares/seqGraph/build/matching", help="the matching software, default is in system env.", metavar="\b")
    optional.add_argument("-h", "--help", action="help")

    options = parser.parse_args()

    if len(sys.argv) == 1:
        print (f"see python {sys.argv[0]} -h")
    else:

        min_q = 20  # read quality cutoff
        min_dp_ratio = 0.1 #0.1 #0.3
        min_map_ratio = 0.3 # for each read

        seg_ref = options.seg_ref
        bam_prefix = options.o + "/" + options.s + ".seg"
        split_graph = options.o + "/" + options.s + ".split.graph"  # graph constructed by only split reads
        pair_graph = options.o + "/" + options.s + ".pair.graph"  # graph constructed by paired-end reads
        bam_name = bam_prefix + ".bam"
        depth_file = bam_prefix + ".bam.depth"
        gene_graph_flag = False
        if options.g == 1 and options.ref_genome != '':
            gene_graph_flag = True
        
        alignment(options.align_tool, options.seg_ref, options.fq1, options.fq2, bam_prefix, options.t)
        os.system(f"samtools depth -aa {bam_name} >{depth_file}")

        contig_len_dict = rename_contigs(seg_ref)
        print ("start extract edges...")
        chrom_copy, median_depth, chrom_depth = cal_copy_number()
        print ("median depth:", median_depth)
        min_edge_dp = min_dp_ratio * median_depth
        calCrossReads(bam_name, split_graph, 'split')
        calCrossReads(bam_name, pair_graph, 'pair')

        plot_graph(split_graph, split_graph+".pdf")
        plot_graph(pair_graph, pair_graph+".pdf")

        print ("graphs constructed, stored in %s and %s"%(split_graph, pair_graph))
        if gene_graph_flag:
            print ("get gene graph...")
            get_gene_graph(options.ref_genome, options.seg_ref, options.o, options.s, options.t)


        print ("start matching...")
        run_match(gene_graph_flag)
        print ("assembly result is in %s"%(options.o + "/" + options.s + ".contigs.fasta"))
        if options.p == 1:
            print ("refine the contigs...")
            refine_assembly()
            print ("refined assembly result is in %s"%(options.o + "/" + options.s + ".contigs.refine.fasta"))

            