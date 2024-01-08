ref=/mnt/d/breakpoints/assembly/sim/database/Escherichia_coli/NZ_CP049085.2.fasta
truth=/mnt/d/breakpoints/assembly/simulation/assembly_test/sim//fasta/NZ_AP022525.1.fasta
sample=sim
fq1=/mnt/d/breakpoints/assembly/simulation/assembly_test/sim//fastq/NZ_AP022525.1.1.fq
fq2=/mnt/d/breakpoints/assembly/simulation/assembly_test/sim//fastq/NZ_AP022525.1.2.fq

#bash ../scripts/workflow5.sh $sample $fq1 $fq2  output/ $ref $truth
bash ../scripts/workflow6.sh $sample $fq1 $fq2  output/ $ref $truth
