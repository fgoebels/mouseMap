#!/bin/bash
scoress=(poisson Pearson wcc apex Jaccard Euclidiean Herdin MatrixNorms)
scores=(Euclidiean)
for i in "${scores[@]}"
do
submitjob 50 src/CalculateScores.py $1 $i $2.$i.txt 
done
