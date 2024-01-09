# requirments: bwa samtools spades(download binary) lumpyexpress blast freebayes vcftools mummer
# python modules: pysam  pyyaml Biopython scikit-learn pyfaidx matplotlib pandas seaborn networkx
# sudo apt -y install delly
# gnuplot

ref=$1
fq1=$2
fq2=$3
ID=$4
outdir=$5
# true=$6

sample=$outdir/$ID
ins_ref=$outdir/$ID.ins.fa
seg_ref=$outdir/$ID.seg.fa
origin_seg_ref=$outdir/$ID.origin.seg.fa

dir=$(cd `dirname $0`; pwd)
start=$(date +%s)
mkdir $outdir


map_qual=20
threads=15



bwa index $ref
samtools faidx $ref
bwa mem -t $threads -R "@RG\tID:id\tSM:sample\tLB:lib" $ref $fq1 $fq2 \
  | samtools view -bhS -> $sample.unsort.bam
samtools sort -o $sample.init.bam $sample.unsort.bam
samtools index $sample.init.bam
samtools depth $sample.init.bam >$sample.init.depth
# delly call -g $ref -o $sample.delly.sv.bcf $sample.init.bam
# bcftools view $sample.delly.sv.bcf >$sample.delly.sv.vcf
svaba run -t $sample.init.bam -G $ref -a $sample -p $threads
python $dir/extract_short_insertions.py $sample.svaba.indel.vcf $sample.short.INS.fasta $sample.short.INS.pos.txt
# python3 $dir/ref_segment.py $ref $seg_ref $sample.delly.sv.vcf $sample.init.depth $sample.init.bam $sample.short.segs.txt $sample.short.INS.pos.txt
python3 $dir/ref_segment.py $ref $seg_ref $sample.svaba.sv.vcf $sample.init.depth $sample.init.bam $sample.short.segs.txt $sample.short.INS.pos.txt
cat $sample.short.INS.fasta >>$seg_ref

end=$(date +%s)
take=$(( end - start ))
echo one Time taken to map reads is ${take} seconds. # >> ${sample}.log


python3 $dir/extract_unmap.py $sample.init.bam $sample.unmap.bam $map_qual $sample.short.segs.txt
samtools fastq -1 $sample.unmapped.1.fq -2 $sample.unmapped.2.fq -s $sample.unmapped.s.fq -@ 8 $sample.unmap.bam
gzip -f $sample.unmapped.*fq
rm -r $outdir/ass
echo "start spades..."
spades.py --isolate -t $threads -1 $sample.unmapped.1.fq.gz -2 $sample.unmapped.2.fq.gz -s $sample.unmapped.s.fq.gz -o $outdir/ass >$sample.spades.log
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
echo "matching is done..."
/home/wangshuai/softwares/seqGraph/build/matching -b --model 1 -v 1 -g $sample.graph.txt -r $sample.solve.path.txt -c $sample.solve.c.path.txt -m $sample.new.graph.txt --break_c

end=$(date +%s)
take=$(( end - start ))
echo three Time taken to map reads is ${take} seconds. # >> ${sample}.log


python3 $dir/graph2contig.py $seg_ref $sample.solve.path.txt $sample.contigs.fasta
bwa index $sample.contigs.fasta
samtools faidx $sample.contigs.fasta
bwa mem -t $threads -R "@RG\tID:id\tSM:sample\tLB:lib" $sample.contigs.fasta $fq1 $fq2 \
  | samtools view -bhS -> $sample.contigs.unsort.bam
samtools sort -o $sample.contigs.bam $sample.contigs.unsort.bam
samtools index $sample.contigs.bam
!
end=$(date +%s)
take=$(( end - start ))
echo four Time taken to map reads is ${take} seconds. # >> ${sample}.log
freebayes -f $sample.contigs.fasta -p 1 $sample.contigs.bam >$sample.contigs.vcf
# delly call -g $sample.contigs.fasta -o $sample.new.sv.bcf $sample.contigs.bam
# bcftools view $sample.new.sv.bcf >$sample.new.sv.vcf
svaba run -t $sample.contigs.bam -G $sample.contigs.fasta -a $sample.svaba_2_ -p $threads
python $dir/polish_contigs.py $sample.contigs.vcf $sample.svaba_2_.svaba.indel.vcf $sample.svaba_2_.svaba.sv.vcf $sample.contigs.fasta $sample.contigs.final.fasta
end=$(date +%s)
take=$(( end - start ))
echo five Time taken to map reads is ${take} seconds. # >> ${sample}.log
# vcftools --vcf $sample.contigs.vcf --minQ 20 --recode --recode-INFO-all --out $sample.contigs_q20
# bgzip -f $sample.contigs_q20.recode.vcf
# tabix -f $sample.contigs_q20.recode.vcf.gz
# cat $sample.contigs.fasta |bcftools consensus -H 1 $sample.contigs_q20.recode.vcf.gz >$sample.contigs.consensus.fasta



# bwa index $sample.contigs.consensus.fasta
# samtools faidx $sample.contigs.consensus.fasta
# bwa mem -t $threads -R "@RG\tID:id\tSM:sample\tLB:lib" $sample.contigs.consensus.fasta $fq1 $fq2 \
#   | samtools view -bhS -> $sample.contigs.new.unsort.bam
# samtools sort -o $sample.contigs.new.bam $sample.contigs.new.unsort.bam
# samtools index $sample.contigs.new.bam
# svaba run -t $sample.contigs.new.bam -G $sample.contigs.consensus.fasta -a $sample.svaba_new_
# python $dir/split_contigs.py $sample.svaba_new_.svaba.sv.vcf $sample.contigs.consensus.fasta $sample.contigs.consensus.split.fasta



end=$(date +%s)
take=$(( end - start ))
echo six Time taken to map reads is ${take} seconds. # >> ${sample}.log


# blastn -subject $seg_ref -query $seg_ref -out $sample.blast.txt -outfmt 6 #remove the overlap between contig and ref-segments
# python $dir/delete_repeat.py $seg_ref $seg_ref.new $sample.blast.txt
# mv $seg_ref.new $seg_ref