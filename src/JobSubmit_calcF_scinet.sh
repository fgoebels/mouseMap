#!/bin/bash

#PBS -l nodes=1:ppn=8,walltime=00:20:00
#PBS -N calcF-Candid
#PBS -q debug

module load intel/15.0.2 R/3.1.1 python/2.7.8

cd /home/g/gbader/fgoebels/mouseMap

echo $INPUT
echo $OUTPUT

src/Calcscores.sh Yeast/data/Yeast_240.txt src/TCSS/gene_ontology.obo.txt Yeast/data/gene_association.tab $INPUT $OUTPUT
