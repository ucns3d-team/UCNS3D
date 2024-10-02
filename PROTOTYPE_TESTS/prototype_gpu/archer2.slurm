#!/bin/bash


#SBATCH --job-name=UCNS3D_PROTOTYPE_TESTING
#SBATCH --gpus=4
#SBATCH --time=00:60:00
#SBATCH --nodes=1
#SBATCH --exclusive

#SBATCH --partition=gpu
#SBATCH --qos=gpu-exc
#SBATCH --cpus-per-task=1
#SBATCH --gpu-bind=single:1

#SBATCH --account=ecseak03

# Setup the job environment (this module needs to be loaded before any other modules)
module load rocm
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

export LD_LIBRARY_PATH=/work/ecseak03/ecseak03/shared/cce_build/local/lib:$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH

export OMP_NUM_THREADS=1
export OMP_PLACES=cores

#export CRAY_ACC_USE_UNIFIED_MEM=1
#export HSA_XNACK=1


srun --ntasks=4 --cpus-per-task=1 --hint=nomultithread --distribution=block:block ./prototype.x




