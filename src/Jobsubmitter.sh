#!/bin/bash
while read l; do
	for NUM in `seq 1 1 120`; do
		qsub src/JobSubmit_calcF_scinet.sh -v INPUT="src/GetEntropyPerFrac.py $l data/WORM/Ce_goldstandart.txt $NUM /scratch/g/gbader/fgoebels/window_out/Worm_1D"
	done
done < $HOME/mouseMap/data/WORM/Allfracs.txt
