ID=$1
fq1=$2
fq2=$3
outdir=$4
sample=$outdir/$ID/$ID
true=$5

mkdir $outdir/$ID

spades.py -t 5 --pe1-1 $fq1 --pe1-2 $fq2 -o $outdir/$ID

mummer -mumreference -l 1000 -b -c $true $sample.contigs.consensus.fasta >$sample.mums
mummerplot -png -p $sample $sample.mums
python $dir/measure.py $sample.mums $true $sample
