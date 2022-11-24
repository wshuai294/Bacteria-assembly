import pysam
import sys

inbam = sys.argv[1]
outbam = sys.argv[2]
q = int(sys.argv[3])
mapbam = sys.argv[4]


samfile = pysam.AlignmentFile(inbam, "rb")
filter_file = pysam.AlignmentFile(outbam, "wb", template=samfile)
filter_origin_file = pysam.AlignmentFile(mapbam, "wb", template=samfile)
for read in samfile.fetch(until_eof=True):
    # if read.is_unmapped or read.mapping_quality < q:
    # because low quality may only indicate multiple alignments
    # so set quality cutoff may bring some multi-aligned reads
    if read.is_unmapped : 
        filter_file.write(read)
    else:
        filter_origin_file.write(read)

 
filter_origin_file.close()
filter_file.close()
samfile.close()
