#!/bin/bash

# Slurm job options (job-name, compute nodes, job time)
#SBATCH --job-name=UCNS3D_PROTOTYPE_TESTING
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --tasks-per-node=2
#SBATCH --cpus-per-task=64

# Replace [budget code] below with your project code (e.g. t01)
#SBATCH --account=ecseak03
#SBATCH --partition=standard
#SBATCH --qos=standard

export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK


# Set the number of threads to 16 and specify placement
#   There are 16 OpenMP threads per MPI process
#   We want one thread per physical core
export OMP_NUM_THREADS=64
export OMP_PLACES=cores

# Launch the parallel job
#   Using 32 MPI processes
#   8 MPI processes per node
#   16 OpenMP threads per MPI process
#   Additional srun options to pin one thread per physical core



srun --hint=nomultithread --distribution=block:block ./prototype.x




