ref=$1
fq1=$2
fq2=$3
ID=$4
outdir=$5
sample=$outdir/$ID

dir=$(cd `dirname $0`; pwd)
start=$(date +%s)
:<<!

###########################detect breakpoints###########################
bwa index $ref
samtools faidx $ref
bwa mem -M -t 5 -R "@RG\tID:id\tSM:sample\tLB:lib" $ref $fq1 $fq2 \
  | samtools view -bhS -> $sample.unsort.bam
samtools sort -o $sample.bam $sample.unsort.bam
rm $sample.unsort.bam

samtools view -h $sample.bam \
  | python3 $dir/extractSplitReads_BwaMem.py -i stdin \
  | samtools view -Sb > $sample.unsort.splitters.bam
samtools sort -o $sample.splitters.bam $sample.unsort.splitters.bam
samtools view -q 20 -b $sample.bam > $sample.unique.bam
samtools index $sample.unique.bam
samtools index $sample.splitters.bam
python3 $dir/get_raw_bkp.py -n 0 -u $sample.unique.bam -o $sample.raw.txt
python3 $dir/accurate_bkp.py -n 0 -g $ref -u $sample.unique.bam \
-s $sample.splitters.bam -a $sample.raw.txt -o $sample.acc.txt
###########################detect breakpoints###########################
!

##########################assemble unmapped#################################
samtools view -f 4 -bS $sample.bam > $sample.unmap.bam
samtools fastq -1 $sample.unmapped.1.fq -2 $sample.unmapped.2.fq -s $sample.unmapped.s.fq -@ 8 $sample.unmap.bam
gzip -f $sample.unmapped.*fq 
/home/wangshuai/miniconda3/bin/megahit --tmp-dir /tmp -1 $sample.unmapped.1.fq.gz -2 $sample.unmapped.2.fq.gz -r $sample.unmapped.s.fq.gz -t 8 -o $outdir/ass 

seg_ref=$ref.seg.fa
python3 $dir/ref_segment.py $ref $seg_ref $sample.acc.txt
cat $outdir/ass/final.contigs.fa >>$seg_ref

###########################bi matching###########################
bwa index $seg_ref
samtools faidx $seg_ref
bwa mem -M -t 5 -R "@RG\tID:id\tSM:sample\tLB:lib" $seg_ref $fq1 $fq2 \
  | samtools view -bhS -> $sample.seg.unsort.bam
samtools sort -o $sample.seg.bam $sample.seg.unsort.bam
samtools index $sample.seg.bam

$dir/../bi_matching/build/overlap_full $sample.seg.bam $sample.out.txt
$dir/../bi_matching/build/bi_solve $sample.out.txt $sample.solve.txt
###########################bi matching###########################

end=$(date +%s)
take=$(( end - start ))
echo Time taken to map reads is ${take} seconds.  >> ${sample}.log