ID=$1
fq1=$2
fq2=$3
outdir=$4
ref=$5
truth=$6

sample=$outdir/$ID
dir=$(cd `dirname $0`; pwd)
threads=10

if [ ! -d "$outdir" ]; then
  mkdir $outdir
fi

split_ref=$sample.split.fasta
seg_ref=$sample.seg.fasta

# :<<!
$dir//src/continuous  $ref $fq1 $sample.read.tab 50
python /mnt/d/breakpoints/assembly/scripts/src/find_continous.py $ref $sample.read.tab $sample.split.fasta 1000
python $dir/filter_contig_by_len.py $sample.split.fasta $sample.split.filter.fasta 1000

split_ref=$sample.split.filter.fasta

$dir/../tools/novoindex $split_ref.ndx $split_ref
samtools faidx $split_ref
$dir/../tools/novoalign -d $split_ref.ndx -f $fq1 $fq2 -F STDFQ -o SAM \
    -o FullNW --mCPU $threads | samtools view -Sb - | samtools sort -  > $sample.split.bam
samtools index $sample.split.bam
samtools view -f 4 -b $sample.split.bam > $sample.unmap.bam
samtools fastq -1 $sample.unmapped.1.fq -2 $sample.unmapped.2.fq -s $sample.unmapped.s.fq $sample.unmap.bam
gzip -f $sample.unmapped.*fq
rm -r $outdir/ass
# !
cat $split_ref>$seg_ref


echo "start spades..."
spades.py --isolate -t $threads -1 $sample.unmapped.1.fq.gz -2 $sample.unmapped.2.fq.gz -s $sample.unmapped.s.fq.gz -o $outdir/ass >$sample.spades.log

if [ -f "$outdir/ass/contigs.fasta" ]; then
    python $dir/filter_assemblies.py $outdir/ass/contigs.fasta $outdir/ass/contigs.filter.fasta 100 10
    contig=$outdir/ass/contigs.filter.fasta
    $dir/../tools/novoindex $contig.ndx $contig
    samtools faidx $contig
    $dir/../tools/novoalign -d $contig.ndx -f $fq1 $fq2 -F STDFQ -o SAM \
        -o FullNW --mCPU $threads | samtools view -Sb - | samtools sort -  > $sample.contig.bam
    samtools index $sample.contig.bam
    python $dir/find_INS_breakpoints.py $sample.contig.bam $contig $outdir/ass/contigs.splitted.fasta
    cat $outdir/ass/contigs.splitted.fasta >>$seg_ref
fi


# cat $seg_ref >$sample.contigs.fasta
# minimap2 -c $seg_ref  $seg_ref > $sample.minimap.paf
# blastn -subject $seg_ref -query $seg_ref -out $sample.blast.txt -outfmt 6 #remove the overlap between contig and ref-segments
# python $dir/delete_repeat.py $seg_ref $seg_ref.uniq.fa $sample.blast.txt
# seg_ref=$seg_ref.uniq.fa
# cat $seg_ref >$sample.contigs.fasta


$dir/../tools/novoindex $seg_ref.ndx $seg_ref
samtools faidx $seg_ref
$dir/../tools/novoalign -d $seg_ref.ndx -f $fq1 $fq2 -F STDFQ -o SAM \
    -o FullNW --mCPU $threads | samtools view -Sb - | samtools sort -  > $sample.seg.bam

# bwa index $seg_ref
# samtools faidx $seg_ref
# bwa mem -a -t $threads -R "@RG\tID:id\tSM:sample\tLB:lib" $seg_ref $fq1 $fq2 | samtools view -Sb - | samtools sort -  > $sample.seg.bam

samtools index $sample.seg.bam
samtools depth -aa $sample.seg.bam >$sample.seg.bam.depth
# rm $sample.seg.unsort.bam
python3 $dir/bam2graph.py $sample.seg.bam $sample.graph.txt $sample.seg.bam.depth $seg_ref
python3 $dir/bam2graph_paired.py $sample.seg.bam $sample.graph_pair.txt $sample.seg.bam.depth $seg_ref
# python3 $dir/plot_graph.py $sample.graph.txt $sample.plot.graph.pdf

###### graph by gene
python3 $dir/gene_distance.py $ref $outdir $seg_ref $ID $threads
sort -k2,2 -k6,6n $sample.gene.graph.txt > $sample.gene.sorted.lh
!
# echo "matching..."
# /home/wangshuai/softwares/seqGraph/build/matching -b --model 1 -v 1 -g $sample.graph.txt -r $sample.solve.path.txt -c $sample.solve.c.path.txt -m $sample.new.graph.txt --break_c --gene_junc $sample.gene.sorted.lh
$dir/../tools/matching -b --model 1 -v 1 -g $sample.graph.txt -r $sample.solve.path.txt -c $sample.solve.c.path.txt -m $sample.new.graph.txt --break_c --gene_junc $sample.gene.sorted.lh
# python3 $dir/graph2contig.py $seg_ref $sample.solve.path.txt $sample.contigs.fasta


# python3 $dir/parse_span.py $sample.graph.txt $sample.graph_pair.txt >$sample.span.junc
# python3 $dir/link_span.py $sample.solve.path.txt $sample.span.junc > $sample.solve.update.path.txt
# python3 $dir/graph2contig.py $seg_ref $sample.solve.update.path.txt $sample.contigs.fasta

#  minimap2 -c $truth $seg_ref > $sample.real.path.paf

cat $seg_ref >$sample.contigs.fasta
rm -r ${sample}_quast
$dir/../tools/quast-5.2.0/quast.py -t 10 $sample.contigs.fasta -r $truth -o ${sample}_quast --strict-NA 
cat ${sample}_quast/report.txt

# rm $sample.contig.bam*
# rm $sample.split.bam*
# rm $sample.unmap.bam*
# rm -r $outdir/ass
# rm $sample.split.fasta
# rm $sample.split.filter.fasta*
# rm ${sample}_ref*
# rm ${sample}_seg*