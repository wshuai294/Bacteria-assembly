ID=$1
fq1=$2
fq2=$3
outdir=$4
sample=$outdir/$ID/$ID
true=$5
dir=$(cd `dirname $0`; pwd)


rm -r $outdir/$ID
spades.py --isolate -t 10 --pe1-1 $fq1 --pe1-2 $fq2 -o $outdir/$ID

# mummer -mumreference -l 1000 -b -c $true $outdir/$ID/contigs.fasta >$sample.mums
# mummerplot -png -Q $outdir/$ID/contigs.fasta -p $sample $sample.mums
# python $dir/measure.py $sample.mums $true $sample $sample.revise.mums


# python3 $dir/measure_blast.py $true $sample $outdir/$ID/contigs.fasta
