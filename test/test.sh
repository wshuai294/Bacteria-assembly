ref=NZ_CP028685.1.fasta
truth=GCF_000008865.2_ASM886v2_genomic.fna
sample=test
fq1=DRR198804_1_subsampled.fastq
fq2=DRR198804_2_subsampled.fastq
gzip -d *gz
#bash ../scripts/workflow5.sh $sample $fq1 $fq2  output/ $ref $truth
python ../scripts/main.py --fq1 DRR198804_1_subsampled.fastq --fq2 DRR198804_2_subsampled.fastq -s test -o output --genome_list genome.list
