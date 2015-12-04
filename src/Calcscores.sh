#!/bin/bash
scores=(Pearson wcc apex Jaccard Euclidiean Herdin MatrixNorms GOSim)
files=()
for i in "${scores[@]}"
do
echo "Calculating score $i"
src/CalculateScores.py $1 $2 $3 $4 $i $5.$i.tmp.txt
cut -f 3- $5.$i.tmp.txt > $5.$i.txt
rm $5.$i.tmp.txt
files+=("$5.$i.txt")
done
paste $4 "${files[@]}" > $5.scores.txt
rm "${files[@]}"
