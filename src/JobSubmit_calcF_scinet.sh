#!/bin/bash

#PBS -l nodes=1:ppn=8,walltime=04:00:00
#PBS -N Worm_window
#PBS -o /scratch/g/gbader/fgoebels/grid.o
#PBS -e /scratch/g/gbader/fgoebels/grid.e

module load intel/15.0.2 R/3.1.1 python/2.7.8

cd /home/g/gbader/fgoebels/mouseMap

eval $INPUT

