#!/bin/bash
#SBATCH --partition=cpsc424

#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --job-name=Serial
#SBATCH --time=1:00:00
#        #SBATCH --exclusive
# ONLY USE THE ABOVE LINE FOR FINAL TIMING RUNS

module load iomkl
pwd
echo $SLURM_JOB_NODELIST
echo $SLURM_NTASKS_PER_NODE
make clean
make serial
time ./serial
time ./serial
time ./serial
