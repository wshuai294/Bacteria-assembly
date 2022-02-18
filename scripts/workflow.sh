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

bwa_score=60 #default 30
map_qual=20

# :<<!

##########################assemble unmapped#################################
bwa index $ref
samtools faidx $ref
bwa mem -T $bwa_score -M -t 5 -R "@RG\tID:id\tSM:sample\tLB:lib" $ref $fq1 $fq2 \
  | samtools view -bhS -> $sample.unsort.bam
samtools sort -o $sample.init.bam $sample.unsort.bam
# samtools view -f 4 -bS $sample.init.bam > $sample.unmap.bam
samtools index $sample.init.bam
python3 $dir/extract_unmap.py $sample.init.bam $sample.unmap.bam $map_qual
samtools fastq -1 $sample.unmapped.1.fq -2 $sample.unmapped.2.fq -s $sample.unmapped.s.fq -@ 8 $sample.unmap.bam
gzip -f $sample.unmapped.*fq 
rm -r $outdir/ass
/mnt/d/breakpoints/tools/SPAdes-3.15.3/bin/spades.py --isolate -t 5 --pe1-1 $sample.unmapped.1.fq.gz\
 --pe1-2 $sample.unmapped.2.fq.gz --pe1-s $sample.unmapped.s.fq.gz -o $outdir/ass >$sample.spades.log
cat $ref $outdir/ass/contigs.fasta >$ins_ref
samtools depth $sample.init.bam >$sample.init.depth
rm $sample.unsort.bam
rm $sample.unmap.bam
# rm $sample.init.bam
# rm $sample.unmapped*

##########################assemble unmapped#################################




# :<<!
##########################lumpy#################################
# Align the data
bwa index $ins_ref
samtools faidx $ins_ref
bwa mem -T $bwa_score -t 5 -R "@RG\tID:id\tSM:sample\tLB:lib" $ins_ref $fq1 $fq2 \
    | samblaster --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20\
    | samtools view -S -b - \
    > $sample.bam
samtools view -b -F 1294 $sample.bam > $sample.discordants.unsorted.bam
samtools view -h $sample.bam \
    | /home/wangshuai/tools/lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin \
    | samtools view -Sb - \
    > $sample.splitters.unsorted.bam
samtools sort $sample.bam -o $sample.sort.bam
samtools sort $sample.discordants.unsorted.bam -o $sample.discordants.bam
samtools sort $sample.splitters.unsorted.bam -o $sample.splitters.bam

samtools view -b -q $map_qual $sample.sort.bam > $sample.sort.unique.bam
samtools depth $sample.sort.unique.bam >$sample.bam.depth
samtools index $sample.sort.unique.bam

lumpyexpress \
    -B $sample.bam \
    -S $sample.splitters.bam \
    -D $sample.discordants.bam \
    -o $sample.sv.vcf
python3 $dir/ref_segment.py $ins_ref $origin_seg_ref $sample.sv.vcf $sample.bam.depth 1 $sample.sort.unique.bam

!

makeblastdb -in $origin_seg_ref -dbtype nucl -out $origin_seg_ref.db -parse_seqids
blastn -query $origin_seg_ref -db $origin_seg_ref.db -outfmt 6 -out $origin_seg_ref.blast.out
python3 $dir/delete_repeat.py $origin_seg_ref $seg_ref $origin_seg_ref.blast.out
rm $ins_ref.*
# rm $sample.bam
rm $sample.splitters*.bam
rm $sample.discordants*.bam
# rm $sample.sort.bam
# rm $sample.bam.depth
##########################lumpy#################################


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


# :<<!
###########################bi matching###########################
seg_num=$(awk '/>/ {a++} END {print a}' $seg_ref)
if [ $seg_num == 1 ]
then
    echo "-----------------there is a single contig.---------------------"
    echo "-----------------there is a single contig.---------------------"
    echo "-----------------there is a single contig.---------------------"
    echo ">contig_1" > $sample.contigs.fasta
    cat $seg_ref |grep -v ">" >> $sample.contigs.fasta
else

    bwa index $seg_ref
    samtools faidx $seg_ref
    bwa mem -M -t 5 -R "@RG\tID:id\tSM:sample\tLB:lib" $seg_ref $fq1 $fq2 \
      | samtools view -bhS -> $sample.seg.unsort.bam
    samtools sort -o $sample.seg.bam $sample.seg.unsort.bam
    samtools index $sample.seg.bam

    echo graph-building...
    $dir/../seqGraph/build/rDistance $sample.seg.bam $sample.raw.graph.txt
    # dp=$(samtools depth $sample.seg.bam  |  awk '{sum+=$3} END { print sum/NR}')
    # python3 $dir/filter_edge.py $sample.raw.graph.txt $sample.graph.txt $dp 0.2
    python3 $dir/bam2graph.py $sample.seg.bam $sample.graph.txt
    # $dir/../seqGraph/build/matching $sample.graph.txt $sample.solve.match.txt $sample.solve.c.path.txt 10 0
    # python3 $dir/bi_matching.py
    python3 $dir/plot_graph.py $sample.graph.txt $sample.plot.graph.pdf
    python3 $dir/graph2contig.py $seg_ref $sample.solve.path.txt $sample.contigs.fasta

    rm $sample.seg.unsort.bam
    # rm $sample.seg.bam*
    # rm $seg_ref.*

fi

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

rm $sample.contigs.fasta.*
rm $sample.contigs.unsort.bam
rm $sample.contigs.bam*
rm $sample.contigs.vcf
###########################consensus sequence###########################


!
mummer -mumreference -l 1000 -b -c $true $sample.contigs.consensus.fasta >$sample.mums

python3 $dir/measure.py $sample.mums $true $sample $sample.revise.mums
mummerplot -png -Q $sample.contigs.consensus.fasta -p $sample $sample.mums

rm $sample.fplot
rm $sample.gp
rm $sample.rplot
rm $sample.contigs_q20.log

makeblastdb -in $true -dbtype nucl -out $true.db -parse_seqids
blastn -query $sample.contigs.consensus.fasta -db $true.db -outfmt 7 -out $sample.map2true.out


end=$(date +%s)
take=$(( end - start ))
echo Time taken to map reads is ${take} seconds.  >> ${sample}.log

!





























































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