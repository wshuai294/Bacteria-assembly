# requirments: bwa samtools spades(download binary) lumpyexpress blast freebayes vcftools mummer
# python modules: pysam  pyyaml Biopython scikit-learn pyfaidx matplotlib pandas seaborn networkx
# sudo apt -y install delly
# gnuplot

ref=$1
fq1=$2
fq2=$3
ID=$4
outdir=$5
true=$6

sample=$outdir/$ID
ins_ref=$outdir/$ID.ins.fa
seg_ref=$outdir/$ID.seg.fa
origin_seg_ref=$outdir/$ID.origin.seg.fa

dir=$(cd `dirname $0`; pwd)
start=$(date +%s)
mkdir $outdir

bwa_score=30 #default 30
map_qual=20
threads=15

# :<<!
bowtie2-build $ref $ref
samtools faidx $ref
bowtie2 -p $threads -x $ref -1 $fq1 -2 $fq2 \
  | samtools view -bhS -> $sample.unsort.bam
samtools sort -o $sample.init.bam $sample.unsort.bam
samtools index $sample.init.bam
samtools depth $sample.init.bam >$sample.init.depth
delly call -g $ref -o $sample.delly.sv.bcf $sample.init.bam
bcftools view $sample.delly.sv.bcf >$sample.delly.sv.vcf
python3 $dir/ref_segment.py $ref $seg_ref $sample.delly.sv.vcf $sample.init.depth $sample.init.bam $sample.short.segs.txt


python3 $dir/extract_unmap.py $sample.init.bam $sample.unmap.bam $map_qual $sample.short.segs.txt
samtools fastq -1 $sample.unmapped.1.fq -2 $sample.unmapped.2.fq -s $sample.unmapped.s.fq -@ 8 $sample.unmap.bam
gzip -f $sample.unmapped.*fq 
rm -r $outdir/ass
echo "start spades..."
spades.py  -t $threads -1 $sample.unmapped.1.fq.gz -2 $sample.unmapped.2.fq.gz -s $sample.unmapped.s.fq.gz --isolate -o $outdir/ass >$sample.spades.log
python $dir/filter_assemblies.py $outdir/ass/contigs.fasta $outdir/ass/contigs.filter.fasta 300
cat $outdir/ass/contigs.filter.fasta >>$seg_ref



bowtie2-build $seg_ref $seg_ref
samtools faidx $seg_ref
bowtie2 -p $threads -x $seg_ref -1 $fq1 -2 $fq2 \
  | samtools view -bhS -> $sample.seg.unsort.bam
samtools sort -o $sample.seg.bam $sample.seg.unsort.bam
samtools index $sample.seg.bam
samtools depth -aa $sample.seg.bam >$sample.seg.bam.depth
rm $sample.seg.unsort.bam

!
echo graph-building...
python3 $dir/bam2graph.py $sample.seg.bam $sample.graph.txt $sample.seg.bam.depth
# :<<!
python3 $dir/plot_graph.py $sample.graph.txt $sample.plot.graph.pdf
echo "matching is done..."
/home/wangshuai/softwares/seqGraph/build/matching -b --model 1 -v 1 -g $sample.graph.txt -r $sample.solve.path.txt -c $sample.solve.c.path.txt -m $sample.new.graph.txt --break_c


python3 $dir/graph2contig.py $seg_ref $sample.solve.path.txt $sample.contigs.fasta
bowtie2-build  $sample.contigs.fasta $sample.contigs.fasta
samtools faidx $sample.contigs.fasta
bowtie2 -p $threads -x $sample.contigs.fasta -1 $fq1 -2 $fq2 \
  | samtools view -bhS -> $sample.contigs.unsort.bam
samtools sort -o $sample.contigs.bam $sample.contigs.unsort.bam
samtools index $sample.contigs.bam
freebayes -f $sample.contigs.fasta -p 1 $sample.contigs.bam >$sample.contigs.vcf
vcftools --vcf $sample.contigs.vcf --minQ 20 --recode --recode-INFO-all --out $sample.contigs_q20
bgzip -f $sample.contigs_q20.recode.vcf
tabix -f $sample.contigs_q20.recode.vcf.gz
cat $sample.contigs.fasta |bcftools consensus -H 1 $sample.contigs_q20.recode.vcf.gz >$sample.contigs.consensus.fasta

rm $sample.contigs.fasta.*
rm $sample.contigs.unsort.bam

mummer -maxmatch -l 1000 -b -c $true $sample.contigs.consensus.fasta >$sample.mums
python3 $dir/measure_blast.py $true $sample $sample.contigs.consensus.fasta

end=$(date +%s)
take=$(( end - start ))
echo Time taken to map reads is ${take} seconds.  >> ${sample}.log

!

