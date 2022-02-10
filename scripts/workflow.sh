ref=$1
fq1=$2
fq2=$3
ID=$4
outdir=$5
true=$6

sample=$outdir/$ID
ins_ref=$outdir/$ID.ins.fa
seg_ref=$outdir/$ID.seg.fa

dir=$(cd `dirname $0`; pwd)
start=$(date +%s)

# :<<!

##########################assemble unmapped#################################
bwa index $ref
samtools faidx $ref
bwa mem -M -t 5 -R "@RG\tID:id\tSM:sample\tLB:lib" $ref $fq1 $fq2 \
  | samtools view -bhS -> $sample.unsort.bam
samtools sort -o $sample.bam $sample.unsort.bam

# samtools index $sample.bam   ###consensus
# freebayes -f $ref -p 1 $sample.bam >$sample.vcf
# vcftools --vcf $sample.vcf --minQ 20 --recode --recode-INFO-all --out $sample.q20
# bgzip -f $sample.q20.recode.vcf
# tabix -f $sample.q20.recode.vcf.gz
# cat $ref |bcftools consensus -H 1 $sample.q20.recode.vcf.gz >$sample.consensus.fasta


samtools view -f 4 -bS $sample.bam > $sample.unmap.bam
samtools fastq -1 $sample.unmapped.1.fq -2 $sample.unmapped.2.fq -s $sample.unmapped.s.fq -@ 8 $sample.unmap.bam
gzip -f $sample.unmapped.*fq 
rm -r $outdir/ass
/home/wangshuai/miniconda3/bin/megahit --tmp-dir /tmp -1 $sample.unmapped.1.fq.gz -2 $sample.unmapped.2.fq.gz \
-r $sample.unmapped.s.fq.gz -t 8 -o $outdir/ass --min-contig-len 500
cat $ref $outdir/ass/final.contigs.fa >$ins_ref
# cat $sample.consensus.fasta $outdir/ass/final.contigs.fa >$ins_ref

rm $sample.unsort.bam
rm $sample.unmap.bam
rm $sample.bam
rm $sample.unmapped*
##########################assemble unmapped#################################
# !



# :<<!
##########################lumpy#################################
# Align the data
bwa index $ins_ref
samtools faidx $ins_ref
bwa mem -t 5 -R "@RG\tID:id\tSM:sample\tLB:lib" $ins_ref $fq1 $fq2 \
    | samblaster --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 \
    | samtools view -S -b - \
    > $sample.bam

# Extract the discordant paired-end alignments.
samtools view -b -F 1294 $sample.bam > $sample.discordants.unsorted.bam

# Extract the split-read alignments
samtools view -h $sample.bam \
    | /home/wangshuai/tools/lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin \
    | samtools view -Sb - \
    > $sample.splitters.unsorted.bam

# Sort both alignments
samtools sort $sample.bam -o $sample.sort.bam
samtools sort $sample.discordants.unsorted.bam -o $sample.discordants.bam
samtools sort $sample.splitters.unsorted.bam -o $sample.splitters.bam
samtools depth $sample.sort.bam >$sample.bam.depth

lumpyexpress \
    -B $sample.bam \
    -S $sample.splitters.bam \
    -D $sample.discordants.bam \
    -o $sample.sv.vcf

python3 $dir/ref_segment.py $ins_ref $seg_ref $sample.sv.vcf $sample.bam.depth 1
rm $ins_ref.*
rm $sample.bam
rm $sample.splitters*.bam
rm $sample.discordants*.bam
rm $sample.sort.bam
rm $sample.bam.depth
##########################lumpy#################################



# :<<!
###########################detect breakpoints###########################
# bwa index $ins_ref
# samtools faidx $ins_ref
# bwa mem -M -t 5 -R "@RG\tID:id\tSM:sample\tLB:lib" $ins_ref $fq1 $fq2 \
#   | samtools view -bhS -> $sample.unsort.bam
# samtools sort -o $sample.bam $sample.unsort.bam
# rm $sample.unsort.bam
# samtools depth $sample.bam >$sample.bam.depth

# samtools view -h $sample.bam \
#   | python3 $dir/extractSplitReads_BwaMem.py -i stdin \
#   | samtools view -Sb > $sample.unsort.splitters.bam
# samtools sort -o $sample.splitters.bam $sample.unsort.splitters.bam
# samtools view -q 20 -b $sample.bam > $sample.unique.bam
# samtools index $sample.unique.bam
# samtools index $sample.splitters.bam

# python3 $dir/get_raw_bkp.py -n 0 -u $sample.unique.bam -o $sample.raw.txt
# python3 $dir/accurate_bkp.py -n 0 -g $ins_ref -u $sample.unique.bam \
# -s $sample.splitters.bam -a $sample.raw.txt -o $sample.acc.txt
# # !
# python3 $dir/ref_segment.py $ins_ref $seg_ref $sample.acc.txt $sample.bam.depth 0
###########################detect breakpoints###########################



# makeblastdb -in $true -dbtype nucl -out $true.db -parse_seqids
# blastn -query $seg_ref -db $true.db -outfmt 7 -out $sample.map2true.out

# :<<!
###########################bi matching###########################
bwa index $seg_ref
samtools faidx $seg_ref
bwa mem -M -t 5 -R "@RG\tID:id\tSM:sample\tLB:lib" $seg_ref $fq1 $fq2 \
  | samtools view -bhS -> $sample.seg.unsort.bam
samtools sort -o $sample.seg.bam $sample.seg.unsort.bam
samtools index $sample.seg.bam

echo graph-building...
$dir/../seqGraph/build/rDistance $sample.seg.bam $sample.graph.txt
$dir/../seqGraph/build/matching $sample.graph.txt $sample.solve.path.txt $sample.solve.c.path.txt 10 0
python3 $dir/graph2contig.py $seg_ref $sample.solve.path.txt $sample.contigs.fasta

rm $sample.seg.unsort.bam
rm $sample.seg.bam*
rm $seg_ref.*
###########################bi matching###########################



###########################consensus sequence###########################
bwa index $sample.contigs.fasta
samtools faidx $sample.contigs.fasta
bwa mem -M -t 5 -R "@RG\tID:id\tSM:sample\tLB:lib" $sample.contigs.fasta $fq1 $fq2 \
  | samtools view -bhS -> $sample.contigs.unsort.bam
samtools sort -o $sample.contigs.bam $sample.contigs.unsort.bam
samtools index $sample.contigs.bam
freebayes -f $sample.contigs.fasta -p 1 $sample.contigs.bam >$sample.contigs.vcf
vcftools --vcf $sample.contigs.vcf --minQ 20 --recode --recode-INFO-all --out $sample.contigs_q20
bgzip -f $sample.contigs_q20.recode.vcf
tabix -f $sample.contigs_q20.recode.vcf.gz
cat $sample.contigs.fasta |bcftools consensus -H 1 $sample.contigs_q20.recode.vcf.gz >$sample.contigs.consensus.fasta

rm $sample.contigs.fasta*
rm $sample.contigs.unsort.bam
rm $sample.contigs.bam*
rm $sample.contigs.vcf
###########################consensus sequence###########################



mummer -maxmatch -l 1000 -b -c $true $sample.contigs.consensus.fasta >$sample.mums
mummerplot -png -p $sample $sample.mums
!
# mummerplot -png -R $true -Q $sample.contigs.consensus.fasta -p $sample $sample.mums
python3 $dir/measure.py $sample.mums $true $sample

rm $sample.fplot
rm $sample.gp
rm $sample.rplot
rm $sample.contigs_q20.log


end=$(date +%s)
take=$(( end - start ))
echo Time taken to map reads is ${take} seconds.  >> ${sample}.log
