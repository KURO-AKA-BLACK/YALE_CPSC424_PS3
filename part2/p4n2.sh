#!/bin/bash
#SBATCH --partition=cpsc424

# set total number of MPI processes
#SBATCH --ntasks=4
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
make task5
make task6
make task7

echo "########## task5 (4 processes and 2 nodes) ##########"
time mpirun -n 4 --bind-to socket --rank-by node:span --report-bindings ./task5
time mpirun -n 4 --bind-to socket --rank-by node:span --report-bindings ./task5
time mpirun -n 4 --bind-to socket --rank-by node:span --report-bindings ./task5

echo "########## task6 (4 processes and 2 nodes) ##########"
time mpirun -n 4 --bind-to socket --rank-by node:span --report-bindings ./task6
time mpirun -n 4 --bind-to socket --rank-by node:span --report-bindings ./task6
time mpirun -n 4 --bind-to socket --rank-by node:span --report-bindings ./task6

echo "########## task7 (4 processes and 2 nodes) ##########"
time mpirun -n 4 --bind-to socket --rank-by node:span --report-bindings ./task7
time mpirun -n 4 --bind-to socket --rank-by node:span --report-bindings ./task7
time mpirun -n 4 --bind-to socket --rank-by node:span --report-bindings ./task7

