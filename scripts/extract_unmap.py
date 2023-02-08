import pysam
import sys

inbam = sys.argv[1]
outbam = sys.argv[2]
q = int(sys.argv[3])
short_segs = sys.argv[4]


def main():

    samfile = pysam.AlignmentFile(inbam, "rb")
    filter_file = pysam.AlignmentFile(outbam, "wb", template=samfile)

    # softclipped_file = open("/home/wangshuai/assembly_result/softclipped.fastq", "w")

    for read in samfile.fetch(until_eof=True):
        # if read.is_unmapped or read.mapping_quality < q:
        # because low quality may only indicate multiple alignments
        # so set quality cutoff may bring some multi-aligned reads
        if read.is_unmapped : 
            filter_file.write(read)

        # if read.cigartuples: # chatGPT
        #     start = 0
        #     for cigar_operation, length in read.cigartuples:
        #         if cigar_operation == 4 and length > 50:
        #             end= start + length
        #             softclipped_seq = read.query_sequence[start:end]
        #             softclipped_qual = read.query_qualities[start:end]
        #             print (read.query_sequence[start:end], read.query_qualities[start:end].tostring().decode())
        #             softclipped_file.write("@" + read.query_name + "\n" + read.query_sequence[start:end] + "\n+\n" + read.query_qualities[start:end].tostring().decode() + "\n")

        #             break
        #         start += length
        # # output the softclipped part of reads

    ########## extract the reads mapped to short segments   
    f = open(short_segs, "r")
    for line in f:
        array = line.strip().split()
        chrom, start, end = array[0], int(array[1]), int(array[2])
        for read in samfile.fetch(chrom, start, end):
            filter_file.write(read)

    filter_file.close()
    samfile.close()


if __name__ == "__main__":
    main()
