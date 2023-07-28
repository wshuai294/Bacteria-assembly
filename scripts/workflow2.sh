# requirments: bwa samtools spades(download binary) lumpyexpress blast freebayes vcftools mummer
# python modules: pysam  pyyaml Biopython scikit-learn pyfaidx matplotlib pandas seaborn networkx
# sudo apt -y install delly
# gnuplot

ref_db=$1
fq1=$2
fq2=$3
ID=$4
outdir=$5
truth=$6

sample=$outdir/$ID
ins_ref=$outdir/$ID.ins.fa
seg_ref=$outdir/$ID.seg.fa
origin_seg_ref=$outdir/$ID.origin.seg.fa

dir=$(cd `dirname $0`; pwd)
start=$(date +%s)
mkdir $outdir


map_qual=20
threads=15

$dir/segmentation $fq1 $fq2 $ref_db 28 10 10 $outdir 30 $ID 0.2
cat $outdir/$ID.map.fasta > $seg_ref

echo "python $dir/classify_unmap_reads.py $fq1 $fq2 $outdir $ID"
python $dir/classify_unmap_reads.py $fq1 $fq2 $outdir $ID
gzip -f $sample.unmapped.*fq

# :<<!
rm -r $outdir/ass
echo "start spades..."
echo "spades.py --isolate -t $threads -1 $sample.unmapped.1.fq.gz -2 $sample.unmapped.2.fq.gz  --isolate -o $outdir/ass >$sample.spades.log"  # -s $sample.unmapped.s.fq.gz
spades.py --isolate -t $threads -1 $sample.unmapped.1.fq.gz -2 $sample.unmapped.2.fq.gz  --isolate -o $outdir/ass >$sample.spades.log  # -s $sample.unmapped.s.fq.gz
python $dir/filter_assemblies.py $outdir/ass/contigs.fasta $outdir/ass/contigs.filter.fasta 100
cat $outdir/ass/contigs.filter.fasta >>$seg_ref

end=$(date +%s)
take=$(( end - start ))
echo two Time taken to map reads is ${take} seconds. # >> ${sample}.log


bwa index $seg_ref
samtools faidx $seg_ref
bwa mem -t $threads -R "@RG\tID:id\tSM:sample\tLB:lib" $seg_ref $fq1 $fq2 \
  | samtools view -bhS -> $sample.seg.unsort.bam
samtools sort -o $sample.seg.bam $sample.seg.unsort.bam
samtools index $sample.seg.bam
samtools depth -aa $sample.seg.bam >$sample.seg.bam.depth
rm $sample.seg.unsort.bam

# !
echo graph-building...
python3 $dir/bam2graph.py $sample.seg.bam $sample.graph.txt $sample.seg.bam.depth
# :<<!
python3 $dir/plot_graph.py $sample.graph.txt $sample.plot.graph.pdf
echo "matching..."
# /home/wangshuai/softwares/seqGraph/build/matching -b --debug --model 1 -v 1 -g $sample.graph.txt -r $sample.solve.path.txt -c $sample.solve.c.path.txt -m $sample.new.graph.txt --break_c
python3 $dir/skip_match_test.py $sample.graph.txt $sample.solve.path.txt
end=$(date +%s)
take=$(( end - start ))
echo three Time taken to map reads is ${take} seconds. # >> ${sample}.log

python3 $dir/graph2contig.py $seg_ref $sample.solve.path.txt $sample.contigs.fasta

##### polish the cotigs
bwa index $sample.contigs.fasta
samtools faidx $sample.contigs.fasta
bwa mem -t $threads -R "@RG\tID:id\tSM:sample\tLB:lib" $sample.contigs.fasta $fq1 $fq2 \
  | samtools view -F 260 -q 20  -bhS -> $sample.contigs.unsort.bam
samtools sort -o $sample.contigs.bam $sample.contigs.unsort.bam
samtools index $sample.contigs.bam
end=$(date +%s)
take=$(( end - start ))
echo four Time taken to alignment to contigs is ${take} seconds. # >> ${sample}.log

# !
# pilon --genome $sample.contigs.fasta --frags $sample.contigs.bam --output $sample.contigs.polish_1 --chunksize 500000
java -Xmx8G -jar /home/wangshuai/softwares/pilon-1.24.jar --genome $sample.contigs.fasta --frags $sample.contigs.bam --output $sample.contigs.polish_1 --fix bases
cp $sample.contigs.polish_1.fasta $sample.contigs.final.fasta


minimap2 -c $truth $sample.contigs.final.fasta > $sample.contigs.final.fasta.paf
minidot $sample.contigs.final.fasta.paf > $sample.contigs.final.fasta.eps && epstopdf $sample.contigs.final.fasta.eps -o $sample.contigs.final.fasta.pdf

minimap2 -c $truth $ref > $sample.paf
minidot $sample.paf > $sample.eps && epstopdf $sample.eps -o $sample.truth.ref.pdf









# bwa index $sample.contigs.fasta
# samtools faidx $sample.contigs.fasta
# # samtools dict -o $sample.contigs.dict $sample.contigs.fasta
# bwa mem -t $threads -R "@RG\tID:id\tSM:sample\tLB:lib" $sample.contigs.fasta $fq1 $fq2 \
#   | samtools view -bhS -> $sample.contigs.unsort.bam
# samtools sort -o $sample.contigs.bam $sample.contigs.unsort.bam
# samtools index $sample.contigs.bam

# end=$(date +%s)
# take=$(( end - start ))
# echo four Time taken to map reads is ${take} seconds. # >> ${sample}.log
# # !
# freebayes -f $sample.contigs.fasta -p 1 $sample.contigs.bam >$sample.contigs.vcf
# # /home/wangshuai/softwares/gatk-4.3.0.0/gatk HaplotypeCaller -ploidy 1 -I $sample.contigs.bam -O $sample.contigs.vcf -R $sample.contigs.fasta 
# end=$(date +%s)
# take=$(( end - start ))
# echo five Time taken to map reads is ${take} seconds. # >> ${sample}.log
# vcftools --vcf $sample.contigs.vcf --minQ 20 --recode --recode-INFO-all --out $sample.contigs_q20
# bgzip -f $sample.contigs_q20.recode.vcf
# tabix -f $sample.contigs_q20.recode.vcf.gz
# cat $sample.contigs.fasta |bcftools consensus -H 1 $sample.contigs_q20.recode.vcf.gz >$sample.contigs.consensus.fasta

# end=$(date +%s)
# take=$(( end - start ))
# echo six Time taken to map reads is ${take} seconds. # >> ${sample}.log
!