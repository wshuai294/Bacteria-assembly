#!/bin/bash

# requirments: bwa samtools spades(download binary) lumpyexpress blast freebayes vcftools mummer svaba vcflib GNU-parallele
# python modules: pysam  pyyaml Biopython scikit-learn pyfaidx matplotlib pandas seaborn networkx


ref_list=$1
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
threads=10


if [ ! -f "$sample.selected.ref.txt" ]; then
    $dir/select_ref $fq1 $fq2 $ref_list 26 10 $threads $sample.match_rate.csv 5
    highest=$(awk -F',' 'BEGIN { max = 0 } 
            { if($2>max) { max=$2; val=$1; second_col=$2 } } 
            END { print val "," second_col }' $sample.match_rate.csv)
    # highest=$(awk -F',' 'BEGIN { max = 0 } 
    #          { if($4>max) { max=$4; val=$1; second_col=$4 } } 
    #          END { print val "," second_col }' $sample.match_rate.csv)
    ref=$(echo $highest | awk -F',' '{print $1}')
    echo "selected ref is $highest."
    echo "$ref" >$sample.selected.ref.txt
else
ref=$(cat $sample.selected.ref.txt)
echo "selected ref is $ref."
fi

# :<<!
bwa index $ref
samtools faidx $ref
bwa mem -t $threads -R "@RG\tID:id\tSM:sample\tLB:lib" $ref $fq1 $fq2 \
  | samtools view -bhS -> $sample.unsort.bam
samtools sort -o $sample.init.bam $sample.unsort.bam
samtools index $sample.init.bam
samtools depth $sample.init.bam >$sample.init.depth
svaba run -t $sample.init.bam -G $ref -a $sample -p $threads
python3 $dir/extract_short_insertions.py $sample.svaba.indel.vcf $sample.short.INS.fasta $sample.short.INS.pos.txt
python3 $dir/ref_segment.py $ref $seg_ref $sample.svaba.sv.vcf $sample.init.depth $sample.init.bam $sample.short.segs.txt $sample.short.INS.pos.txt $sample
cat $sample.short.INS.fasta >>$seg_ref

end=$(date +%s)
take=$(( end - start ))
echo one Time taken to segmentation is ${take} seconds. # >> ${sample}.log


python3 $dir/extract_unmap.py $sample.init.bam $sample.unmap.bam $map_qual $sample.short.segs.txt
samtools sort -o $sample.unmap.sort.bam $sample.unmap.bam
samtools fastq -1 $sample.unmapped.1.fq.gz -2 $sample.unmapped.2.fq.gz -s $sample.unmapped.s.fq.gz --threads $threads $sample.unmap.sort.bam
if [ -d "$outdir/ass" ]; then
  rm -r $outdir/ass
fi
echo "start spades..."
spades.py --isolate -t $threads -1 $sample.unmapped.1.fq.gz -2 $sample.unmapped.2.fq.gz -s $sample.unmapped.s.fq.gz -o $outdir/ass --cov-cutoff 'auto' >$sample.spades.log
echo "spades.py --isolate -t $threads --12 $sample.unmapped.s.fq.gz -o $outdir/ass --cov-cutoff 'auto' >$sample.spades.log"

python $dir/filter_assemblies.py $outdir/ass/contigs.fasta $outdir/ass/contigs.filter.fasta 1000 # original 100
if [ -f "$outdir/ass/contigs.filter.fasta" ]; then
  cat $outdir/ass/contigs.filter.fasta >>$seg_ref
fi

# blastn -subject $seg_ref -query $seg_ref -out $sample.blast.txt -outfmt 6 #remove the overlap between contig and ref-segments
# python $dir/delete_repeat.py $seg_ref $seg_ref.uniq.fa $sample.blast.txt
# seg_ref=$seg_ref.uniq.fa
# :<<!
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

# !
echo graph-building...
python3 $dir/bam2graph.py $sample.seg.bam $sample.graph.txt $sample.seg.bam.depth
python3 $dir/plot_graph.py $sample.graph.txt $sample.plot.graph.pdf
# python3 $dir/skip_match_test.py $sample.graph.txt $sample.solve.path.txt
/home/wangshuai/softwares/seqGraph/build/matching -b --debug --model 1 -v 1 -g $sample.graph.txt -r $sample.solve.path.txt -c $sample.solve.c.path.txt -m $sample.new.graph.txt --break_c
python3 $dir/graph2contig.py $seg_ref $sample.solve.path.txt $sample.contigs.fasta
echo "matching is done..."
end=$(date +%s)
take=$(( end - start ))
echo three Time taken to matching is ${take} seconds. # >> ${sample}.log


##### polish the cotigs
bwa index $sample.contigs.fasta
samtools faidx $sample.contigs.fasta
bwa mem -t $threads -R "@RG\tID:id\tSM:sample\tLB:lib" $sample.contigs.fasta $fq1 $fq2 \
  | samtools view -F 260 -q 30  -bhS -> $sample.contigs.unsort.bam # -L $sample.solve.path.txt.ref.bed
samtools sort -o $sample.contigs.bam $sample.contigs.unsort.bam
samtools index $sample.contigs.bam
end=$(date +%s)
take=$(( end - start ))
echo four Time taken to alignment to contigs is ${take} seconds. # >> ${sample}.log
!

# pilon --genome $sample.contigs.fasta --frags $sample.contigs.bam --output $sample.contigs.polish_1 --chunksize 500000
java -Xmx8G -jar /home/wangshuai/softwares/pilon-1.24.jar --genome $sample.contigs.fasta --frags $sample.contigs.bam --output $sample.contigs.polish_1 --chunksize 1000000
cp $sample.contigs.polish_1.fasta $sample.contigs.final.fasta
minimap2 -c $truth $sample.contigs.final.fasta > $sample.contigs.final.fasta.paf
minidot $sample.contigs.final.fasta.paf > $sample.contigs.final.fasta.eps && epstopdf $sample.contigs.final.fasta.eps -o $sample.contigs.final.fasta.pdf

minimap2 -c $truth $ref > $sample.paf
minidot $sample.paf > $sample.eps && epstopdf $sample.eps -o $sample.truth.ref.pdf

# python3 $dir/resolve_misassembly.py -i $sample.contigs.polish_1.fasta -o $sample.contigs.final.fasta -r $fq1


# contig_file=$sample.contigs.polish_1.fasta
# bam_file=$sample.contigs.2.bam
# pileup_file=$sample.contigs.samtools.pileup
# bwa index $contig_file
# samtools faidx $contig_file
# bwa mem -a -t $threads $contig_file $fq1 $fq2 | samtools view -h -q 10 -m 50 -F 4 -b | samtools sort > $bam_file
# samtools mpileup -C 50 -A -f $contig_file $bam_file |  awk '$3 != "N"' > $pileup_file
# metaMIC extract_feature --mlen 500 --bam $bam_file -c $contig_file -o $outdir --pileup $pileup_file -m single
# metaMIC predict --mlen 500 -c $contig_file -o $outdir -m single --st 0.3 --slen 100 --nb 2  --at 0.3
# cp $outdir/metaMIC_corrected_contigs.fa $sample.contigs.final.fasta

# ### freebayes-parallele requires the package: parallele and vcflib, can be installed by conda, 
# ### remember to edit the conda_env/bin/vcffirstheader with python2
# freebayes-parallel <(python $dir/fasta_generate_regions.py $sample.contigs.fasta.fai 100000) $threads -p 1 -f $sample.contigs.fasta $sample.contigs.bam >$sample.contigs.vcf
# # freebayes -f $sample.contigs.fasta -t $sample.solve.path.txt.ref.bed -p 1 $sample.contigs.bam >$sample.contigs.vcf
# bgzip -f $sample.contigs.vcf
# tabix -f $sample.contigs.vcf.gz
# svaba run -t $sample.contigs.bam -G $sample.contigs.fasta -a $sample.svaba_2_ -p $threads
# python $dir/polish_contigs.py $sample.contigs.vcf.gz $sample.svaba_2_.svaba.indel.vcf $sample.svaba_2_.svaba.sv.vcf $sample.contigs.fasta $sample.contigs.final.fasta




end=$(date +%s)
take=$(( end - start ))
echo five Time taken to variants is ${take} seconds. # >> ${sample}.log
# !
echo ""
echo ""
echo "**************Assembly Finishied*****************"
echo ""
echo ""
!
