# requirments: bwa samtools spades(download binary) lumpyexpress blast freebayes vcftools mummer svaba vcflib GNU-parallele
# python modules: pysam  pyyaml Biopython scikit-learn pyfaidx matplotlib pandas seaborn networkx


ref_list=$1
fq1=$2
fq2=$3
ID=$4
outdir=$5

sample=$outdir/$ID
ins_ref=$outdir/$ID.ins.fa
seg_ref=$outdir/$ID.seg.fa
origin_seg_ref=$outdir/$ID.origin.seg.fa

dir=$(cd `dirname $0`; pwd)
start=$(date +%s)
mkdir $outdir


map_qual=20
threads=10

$dir/select_ref $fq1 $fq2 $ref_list 26 10 $threads $sample.match_rate.csv 30
# ref=$(awk -F',' 'BEGIN { max = 0 } 
#          { if($2>max) { max=$2; val=$1 } } 
#          END { print val }' $sample.match_rate.csv)
highest=$(awk -F',' 'BEGIN { max = 0 } 
         { if($2>max) { max=$2; val=$1; second_col=$2 } } 
         END { print val "," second_col }' $sample.match_rate.csv)
ref=$(echo $highest | awk -F',' '{print $1}')
echo "selected ref is $highest."

# :<<!
bwa index $ref
samtools faidx $ref
bwa mem -t $threads -R "@RG\tID:id\tSM:sample\tLB:lib" $ref $fq1 $fq2 \
  | samtools view -bhS -> $sample.unsort.bam
samtools sort -o $sample.init.bam $sample.unsort.bam
samtools index $sample.init.bam
samtools depth $sample.init.bam >$sample.init.depth
svaba run -t $sample.init.bam -G $ref -a $sample -p $threads
python $dir/extract_short_insertions.py $sample.svaba.indel.vcf $sample.short.INS.fasta $sample.short.INS.pos.txt
python3 $dir/ref_segment.py $ref $seg_ref $sample.svaba.sv.vcf $sample.init.depth $sample.init.bam $sample.short.segs.txt $sample.short.INS.pos.txt
cat $sample.short.INS.fasta >>$seg_ref

end=$(date +%s)
take=$(( end - start ))
echo one Time taken to segmentation is ${take} seconds. # >> ${sample}.log


python3 $dir/extract_unmap.py $sample.init.bam $sample.unmap.bam $map_qual $sample.short.segs.txt
samtools fastq -1 $sample.unmapped.1.fq -2 $sample.unmapped.2.fq -s $sample.unmapped.s.fq -@ 8 $sample.unmap.bam
gzip -f $sample.unmapped.*fq
if [ -d "$outdir/ass" ]; then
  rm -r $outdir/ass
fi
echo "start spades..."
spades.py --isolate -t $threads -1 $sample.unmapped.1.fq.gz -2 $sample.unmapped.2.fq.gz -s $sample.unmapped.s.fq.gz -o $outdir/ass >$sample.spades.log
python $dir/filter_assemblies.py $outdir/ass/contigs.fasta $outdir/ass/contigs.filter.fasta 100 # original 100
if [ -f "$outdir/ass/contigs.filter.fasta" ]; then
  cat $outdir/ass/contigs.filter.fasta >>$seg_ref
fi

end=$(date +%s)
take=$(( end - start ))
echo two Time taken to unmap-reads assembly is ${take} seconds. # >> ${sample}.log


bwa index $seg_ref
samtools faidx $seg_ref
bwa mem -t $threads -R "@RG\tID:id\tSM:sample\tLB:lib" $seg_ref $fq1 $fq2 \
  | samtools view -bhS -> $sample.seg.unsort.bam
samtools sort -o $sample.seg.bam $sample.seg.unsort.bam
samtools index $sample.seg.bam
samtools depth -aa $sample.seg.bam >$sample.seg.bam.depth
rm $sample.seg.unsort.bam


echo graph-building...
python3 $dir/bam2graph.py $sample.seg.bam $sample.graph.txt $sample.seg.bam.depth
python3 $dir/plot_graph.py $sample.graph.txt $sample.plot.graph.pdf
echo "matching is done..."
/home/wangshuai/softwares/seqGraph/build/matching -b --model 1 -v 1 -g $sample.graph.txt -r $sample.solve.path.txt -c $sample.solve.c.path.txt -m $sample.new.graph.txt --break_c

end=$(date +%s)
take=$(( end - start ))
echo three Time taken to matching is ${take} seconds. # >> ${sample}.log


python3 $dir/graph2contig.py $seg_ref $sample.solve.path.txt $sample.contigs.fasta
bwa index $sample.contigs.fasta
samtools faidx $sample.contigs.fasta
bwa mem -t $threads -R "@RG\tID:id\tSM:sample\tLB:lib" $sample.contigs.fasta $fq1 $fq2 \
  | samtools view -bhS -> $sample.contigs.unsort.bam
samtools sort -o $sample.contigs.bam $sample.contigs.unsort.bam
samtools index $sample.contigs.bam
end=$(date +%s)
take=$(( end - start ))
echo four Time taken to alignment to contigs is ${take} seconds. # >> ${sample}.log


### freebayes-parallele requires the package: parallele and vcflib, can be installed by conda, 
### remember to edit the conda_env/bin/vcffirstheader with python2
freebayes-parallel <(python $dir/fasta_generate_regions.py $sample.contigs.fasta.fai 100000) $threads -p 1 -f $sample.contigs.fasta $sample.contigs.bam >$sample.contigs.vcf
# freebayes -f $sample.contigs.fasta -p 1 $sample.contigs.bam >$sample.contigs.vcf
svaba run -t $sample.contigs.bam -G $sample.contigs.fasta -a $sample.svaba_2_ -p $threads
python $dir/polish_contigs.py $sample.contigs.vcf $sample.svaba_2_.svaba.indel.vcf $sample.svaba_2_.svaba.sv.vcf $sample.contigs.fasta $sample.contigs.final.fasta
end=$(date +%s)
take=$(( end - start ))
echo five Time taken to variants is ${take} seconds. # >> ${sample}.log
# !
echo ""
echo ""
echo "**************Assembly Finishied*****************"
echo ""
echo ""
