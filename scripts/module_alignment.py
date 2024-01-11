import os
import sys



def alignment(alignment_tool, ref, fq1, fq2, prefix, threads):

    header = f"""
        dir={sys.path[0]}
        fq1={fq1}
        fq2={fq2}
        ref={ref}
        threads={threads}
        prefix={prefix}
        
        echo map reads to $ref...

    """

    if alignment_tool.lower() == "novoalign":
        command = f"""

        $dir/../tools/novoindex $ref.ndx $ref
        samtools faidx $ref
        $dir/../tools/novoalign -d $ref.ndx -f $fq1 $fq2 -F STDFQ -o SAM \
            -o FullNW --mCPU $threads | samtools view -Sb - | samtools sort -  > $prefix.bam
        samtools index $prefix.bam
        """
    elif alignment_tool.lower() == "bwa":
        # print ("bwa")
        command = f"""

        bwa index $ref
        samtools faidx $ref

        bwa mem -t $threads -R "@RG\\tID:id\\tSM:sample\\tLB:lib" $ref $fq1 $fq2 \
        | samtools view -Sb - | samtools sort -  > $prefix.bam
        samtools index $prefix.bam
        """
        
    else:
        print ("please select alignment tool")
        sys.exit()

    os.system(header + command)