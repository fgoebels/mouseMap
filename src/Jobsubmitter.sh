#!/bin/bash
for NUM in `seq 1 1 240`;
do
submitjob 10 src/CalculateCoElutionScores.py 2D_Worm/Contrast_Celegans_240.sel.sorted.mapped.txt 2D_Worm/Ce_goldstandard.txt $NUM 2D_Worm/$NUM.test.out
done
