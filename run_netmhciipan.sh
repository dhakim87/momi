#!/bin/bash -e

#PBS -N NetMHCIIPan
#PBS -l walltime=1000:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=1gb
#PBS -V
#PBS -j oe
#PBS -d .
#PBS -t 0-999
#PBS -M djhakim@eng.ucsd.edu

# to run individual jobs, comma-separate the indices
# -t 1,2
# to run all jobs:  
# -t 1-32%16

source activate imsms
NETMHCPATH = "/home/djhakim/netMHCIIpan/netMHCIIpan-4.0/netMHCIIpan"
python run_netmhciipan.py ${PBS_ARRAYID} 1000 $NETMHCPATH

