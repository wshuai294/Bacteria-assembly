import pysam
import sys

inbam = sys.argv[1]
outbam = sys.argv[2]
q = int(sys.argv[3])


samfile = pysam.AlignmentFile(inbam, "rb")
filter_file = pysam.AlignmentFile(outbam, "wb", template=samfile)
for read in samfile.fetch(until_eof=True):
    # if read.flag == 4 or read.mapping_quality < q:
    if read.is_unmapped or read.mapping_quality < q:
        filter_file.write(read)

 

filter_file.close()
samfile.close()
