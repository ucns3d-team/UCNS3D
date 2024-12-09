#!/bin/bash
#
#SBATCH --job-name=ucns3d_prototype
#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --gres=gpu:4
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --cpus-per-task=10
#SBATCH --time=04:00:00
#SBATCH --gpu-bind=single:1
#SBATCH --output=prototypetest_Ne_2e6.out

# Load the required modules
module load openmpi/4.1.6-cuda-12.4-nvfortran
module load nvidia/nvhpc-nompi
module load gcc
#module load ucx/1.16.0
module load intel-20.4/cmkl

export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

export NV_ACC_CUDA_HEAPSIZE=2GB
#export NV_ACC_DEBUG=1
export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export NVCOMPILER_TERM=trace
export OMPI_MCA_BTL=ob1,tcp


srun --ntasks=4 --cpus-per-task=10 ./prototype_gpu.x 2000000
