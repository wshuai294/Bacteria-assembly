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

# :<<!
bwa index $ref
samtools faidx $ref
bwa mem -T $bwa_score -M -t 5 -R "@RG\tID:id\tSM:sample\tLB:lib" $ref $fq1 $fq2 \
  | samtools view -bhS -> $sample.unsort.bam
samtools sort -o $sample.init.bam $sample.unsort.bam
samtools index $sample.init.bam
python3 $dir/extract_unmap.py $sample.init.bam $sample.unmap.bam $map_qual $sample.init.map.bam
samtools fastq -1 $sample.unmapped.1.fq -2 $sample.unmapped.2.fq -s $sample.unmapped.s.fq -@ 8 $sample.unmap.bam
gzip -f $sample.unmapped.*fq 
rm -r $outdir/ass
echo "start spades..."
spades.py  -t 5 -1 $sample.unmapped.1.fq.gz -2 $sample.unmapped.2.fq.gz -s $sample.unmapped.s.fq.gz --isolate -o $outdir/ass >$sample.spades.log

samtools index $sample.init.map.bam
samtools depth $sample.init.map.bam >$sample.init.map.depth
rm $sample.unsort.bam
rm $sample.unmap.bam

delly call -g $ref -o $sample.delly.sv.bcf $sample.init.map.bam
bcftools view $sample.delly.sv.bcf >$sample.delly.sv.vcf

python3 $dir/ref_segment.py $ref $origin_seg_ref $sample.delly.sv.vcf $sample.init.map.depth 1 $sample.init.map.bam
cat $outdir/ass/contigs.fasta >>$origin_seg_ref

makeblastdb -in $origin_seg_ref -dbtype nucl -out $origin_seg_ref.db -parse_seqids
blastn -query $origin_seg_ref -db $origin_seg_ref.db -outfmt 6 -out $origin_seg_ref.blast.out
python3 $dir/delete_repeat.py $origin_seg_ref $seg_ref $origin_seg_ref.blast.out
# rm $ins_ref.*


bwa index $seg_ref
samtools faidx $seg_ref
bwa mem -M -t 5 -R "@RG\tID:id\tSM:sample\tLB:lib" $seg_ref $fq1 $fq2 \
  | samtools view -bhS -> $sample.seg.unsort.bam
samtools sort -o $sample.seg.bam $sample.seg.unsort.bam
samtools index $sample.seg.bam
samtools depth -aa $sample.seg.bam >$sample.seg.bam.depth

echo graph-building...
python3 $dir/bam2graph.py $sample.seg.bam $sample.graph.txt $sample.seg.bam.depth
python3 $dir/plot_graph.py $sample.graph.txt $sample.plot.graph.pdf

/home/wangshuai/softwares/seqGraph/build/matching -b --model 1 -v 1 -g $sample.graph.txt -r $sample.solve.path.txt -c $sample.solve.c.path.txt -m $sample.new.graph.txt --break_c
# $dir/../seqGraph/build/matching --model 1 -g $sample.graph.txt -r $sample.solve.path.txt -c $sample.solve.c.path.txt -i 10 -v 0 -m $sample.new.graph.txt --break_c
echo "matching is done."

python3 $dir/graph2contig.py $seg_ref $sample.solve.path.txt $sample.contigs.fasta
rm $sample.seg.unsort.bam


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

rm $sample.contigs.fasta.*
rm $sample.contigs.unsort.bam

mummer -maxmatch -l 1000 -b -c $true $sample.contigs.consensus.fasta >$sample.mums
# echo "mummer -maxmatch -l 1000 -b -c $true $sample.contigs.consensus.fasta >$sample.mums"

# python3 $dir/measure.py $sample.mums $true $sample $sample.revise.mums
# # :<<!
# mummerplot -png -Q $sample.contigs.consensus.fasta -p $sample $sample.mums
# echo "mummerplot -png -Q $sample.contigs.consensus.fasta -p $sample $sample.mums"

# rm $sample.fplot
# rm $sample.gp
# rm $sample.rplot
# rm $sample.contigs_q20.log
!
makeblastdb -in $true -dbtype nucl -out $true.db -parse_seqids
blastn -query $sample.contigs.consensus.fasta -db $true.db -outfmt 7 -out $sample.map2true.out
python3 $dir/measure_blast.py $sample.map2true.out $true $sample 

end=$(date +%s)
take=$(( end - start ))
echo Time taken to map reads is ${take} seconds.  >> ${sample}.log

# !





























































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

# /home/wangshuai/miniconda3/bin/megahit --tmp-dir /tmp -1 $sample.unmapped.1.fq.gz -2 $sample.unmapped.2.fq.gz \
# -r $sample.unmapped.s.fq.gz -t 8 -o $outdir/ass --min-contig-len 500
# cat $ref $outdir/ass/final.contigs.fa >$ins_ref

# samtools view -h $sample.bam | awk '$1 ~ /@/ ||  $6 ~ /H/{print}' |samtools view -bS - > $sample.hardclip.bam
# samtools merge -f -X $sample.assemble.bam $sample.unmap.bam $sample.hardclip.bam