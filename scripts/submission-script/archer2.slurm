#!/bin/bash

# Slurm job options (job-name, compute nodes, job time)
#SBATCH --job-name=DGUNF5
#SBATCH --time=24:00:00
#SBATCH --nodes=20
#SBATCH --ntasks=160
#SBATCH --tasks-per-node=8
#SBATCH --cpus-per-task=16

# Replace [budget code] below with your project code (e.g. t01)
#SBATCH --account=budget code
#SBATCH --partition=standard
#SBATCH --qos=standard

# Setup the job environment (this module needs to be loaded before any other modules)
module load epcc-job-env

export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH

# Set the number of threads to 16 and specify placement
#   There are 16 OpenMP threads per MPI process
#   We want one thread per physical core
export OMP_NUM_THREADS=16
export OMP_PLACES=cores

# Launch the parallel job
#   Using 32 MPI processes
#   8 MPI processes per node
#   16 OpenMP threads per MPI process
#   Additional srun options to pin one thread per physical core



srun --hint=nomultithread --distribution=block:block ./ucns3d_p




