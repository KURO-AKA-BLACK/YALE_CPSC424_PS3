#!/bin/bash
#SBATCH --partition=cpsc424

# set total number of MPI processes
#SBATCH --ntasks=7
# set number of MPI processes per node
# (number of nodes is calculated by Slurm)
#SBATCH --ntasks-per-node=2
#SBATCH --ntasks-per-socket=1
# set number of cpus per MPI process
#SBATCH --cpus-per-task=1
# set memory per cpu
#SBATCH --mem-per-cpu=5G
#SBATCH --job-name=MPI_RUN
#SBATCH --time=1:00:00

#         #SBATCH --exclusive
# ONLY USE THE ABOVE LINE FOR FINAL TIMING RUNS
module load iomkl
pwd
# echo some environment variables
echo $SLURM_JOB_NODELIST
echo $SLURM_NTASKS_PER_NODE
# Do a clean build
make clean
# My MPI program is task5
make task8

time mpirun -n 7 --bind-to socket --rank-by node:span --report-bindings ./task8
time mpirun -n 7 --bind-to socket --rank-by node:span --report-bindings ./task8
time mpirun -n 7 --bind-to socket --rank-by node:span --report-bindings ./task8

