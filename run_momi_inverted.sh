#!/bin/bash -e

#PBS -N Momi_Inverted
#PBS -l walltime=40:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=1gb
#PBS -V
#PBS -j oe
#PBS -d .
#PBS -t 0-24
#PBS -M djhakim@eng.ucsd.edu

# to run individual jobs, comma-separate the indices
# -t 1,2
# to run all jobs:  
# -t 1-32%16

source activate imsms
python master.py ${PBS_ARRAYID} 25
