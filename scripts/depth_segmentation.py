
depth_file = "/mnt/d/breakpoints/assembly/simulation/test/test.bam.depth"

min_gap = 50

def segmentation(depth_file):
    intervals = []
    pre_chrom = ""
    start = 0
    index = 0

    f = open(depth_file, 'r')
    for line in f:
        array = line.strip().split()
        chrom = array[0]
        pos = int(array[1])

        if pre_chrom == "":
            start = pos
            index = pos
            pre_chrom = chrom
        elif chrom == pre_chrom:
            if pos - index < min_gap:
                index = pos
            else:
                intervals.append([pre_chrom, start, index]) 
                start = pos
                index = pos
        else:
            intervals.append([pre_chrom, start, index]) 
            start = pos
            index = pos
            pre_chrom = chrom

    intervals.append([pre_chrom, start, index]) 
    start = pos
    index = pos
    pre_chrom = chrom

    # print (intervals)
    f.close()


segmentation(depth_file)