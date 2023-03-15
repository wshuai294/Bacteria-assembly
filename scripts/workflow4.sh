#!/bin/bash

# requirments: bwa samtools spades(download binary) lumpyexpress blast freebayes vcftools mummer svaba vcflib GNU-parallele
# python modules: pysam  pyyaml Biopython scikit-learn pyfaidx matplotlib pandas seaborn networkx


ref_list=$1
fq1=$2
fq2=$3
ID=$4
outdir=$5
long=$6

sample=$outdir/$ID
ins_ref=$outdir/$ID.ins.fa
seg_ref=$outdir/$ID.seg.fa
origin_seg_ref=$outdir/$ID.origin.seg.fa

dir=$(cd `dirname $0`; pwd)
start=$(date +%s)
mkdir $outdir


map_qual=20
threads=10

# :<<!
# $dir/select_ref $fq1 $fq2 $ref_list 26 10 $threads $sample.match_rate.csv 30
# highest=$(awk -F',' 'BEGIN { max = 0 } 
#          { if($2>max) { max=$2; val=$1; second_col=$2 } } 
#          END { print val "," second_col }' $sample.match_rate.csv)
# # highest=$(awk -F',' 'BEGIN { max = 0 } 
# #          { if($4>max) { max=$4; val=$1; second_col=$4 } } 
# #          END { print val "," second_col }' $sample.match_rate.csv)
# ref=$(echo $highest | awk -F',' '{print $1}')
# echo "selected ref is $highest."

ref=/mnt/d/breakpoints/assembly/sim/database/ecoli/NZ_CP049085.2.fasta
# ref=test_ref.fasta

pbmm2 align $ref $long $sample.movie1.bam --sort --sample $4 --rg '@RG\tID:movie1'
# pbsv discover $sample.movie1.bam $sample.svsig.gz
# tabix -c '#' -s 3 -b 4 -e 4 $sample.svsig.gz
# pbsv call $ref $sample.svsig.gz $sample.var.vcf
svim alignment $output $sample.movie1.bam $ref --sample $4
python3 $dir/get_consensus_sv.py $ref $sample.var.vcf $sample.contigs.fasta


##### polish the cotigs
bwa index $sample.contigs.fasta
samtools faidx $sample.contigs.fasta
bwa mem -t $threads -R "@RG\tID:id\tSM:sample\tLB:lib" $sample.contigs.fasta $fq1 $fq2 \
  | samtools view -bhS -> $sample.contigs.unsort.bam
samtools sort -o $sample.contigs.bam $sample.contigs.unsort.bam
samtools index $sample.contigs.bam
end=$(date +%s)
take=$(( end - start ))
echo four Time taken to alignment to contigs is ${take} seconds. # >> ${sample}.log
# !
pilon --genome $sample.contigs.fasta --frags $sample.contigs.bam --output $sample.contigs.polish_1 --chunksize 1000000
cp $sample.contigs.polish_1.fasta $sample.contigs.final.fasta


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
