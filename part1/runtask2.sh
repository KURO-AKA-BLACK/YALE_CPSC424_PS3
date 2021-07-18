#!/bin/bash
#SBATCH --partition=cpsc424
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --tasks-per-node=2
#SBATCH --mem-per-cpu=5G
#SBATCH --time=20:00
#SBATCH --job-name=task2
#SBATCH --output=%x-%j.out

# Load necessary module files
module load iomkl
module list

# Print initial working directory
echo " "
echo " "
echo "Working Directory:"
pwd

echo " "
echo " "
echo "Making task2"
make clean
make task2

# Print the node list
echo " "
echo " "
echo "Node List:"
echo $SLURM_NODELIST
echo "ntasks-per-node = " $SLURM_NTASKS_PER_NODE

# Run the program 3 times
echo " "
echo " "
echo "Run 1"
time mpiexec -n 4 task2
echo " "
echo " "
echo "Run 2"
time mpiexec -n 4 task2
echo " "
echo " "
echo "Run 3"
time mpiexec -n 4 task2
