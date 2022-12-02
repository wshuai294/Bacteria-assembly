import pysam
import sys

inbam = sys.argv[1]
outbam = sys.argv[2]
q = int(sys.argv[3])
short_segs = sys.argv[4]


def main():
    samfile = pysam.AlignmentFile(inbam, "rb")
    filter_file = pysam.AlignmentFile(outbam, "wb", template=samfile)
    for read in samfile.fetch(until_eof=True):
        # if read.is_unmapped or read.mapping_quality < q:
        # because low quality may only indicate multiple alignments
        # so set quality cutoff may bring some multi-aligned reads
        if read.is_unmapped : 
            filter_file.write(read)

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
