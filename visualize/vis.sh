minimap2 -c $1 $2 > $1.paf
minidot $1.paf > $1.eps && epstopdf $1.eps -o $3.pdf

rm $1.paf
rm $1.eps
