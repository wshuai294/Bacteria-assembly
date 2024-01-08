ID=$1
fq1=$2
fq2=$3
outdir=$4
ref=$5
truth=$6

sample=$outdir/$ID


/mnt/d/breakpoints/assembly/scripts/src/continuous  $ref $fq1 $sample.read.tab 100

python /mnt/d/breakpoints/assembly/scripts/src/find_continous.py   $ref $sample.read.tab $sample.split.fasta 1000

rm -r ${sample}_quast
~/softwares/quast-5.2.0/quast.py -t 10 $sample.split.fasta -r $truth -o ${sample}_quast --strict-NA 
cat ${sample}_quast/report.txt