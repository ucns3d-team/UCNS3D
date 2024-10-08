#!/bin/bash
#
# UCNS3D : An Open Source, High-Order, Finite-Volume CFD Solver
#
# Main contact:  Panagiotis (Takis) Tsoutsanis (Cranfield University)
#
# Licensed under the GNU General Public License, Version 3.0 (the "License");
# you may not use this file except in compliance with the License.
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "as is" basis,
# without warranties or conditions of any kind, either express or implied.
# see the License for the specific language governing permissions and
# limitations under the License.
#
# --------------------------------------------------------------------
#
# SLURM job submission file for CRAY build on ARCHER2
#

# Slurm job options (job-name, compute nodes, job time)
#SBATCH --job-name=proto_gpu
#SBATCH --output=res_Ne5000.%A
#SBATCH --error=joberr.%A
#SBATCH --time=0:10:0
#SBATCH --account=ecseak03

#SBATCH --partition=gpu
#SBATCH --qos=gpu-exc
#SBATCH --gpus=1
#SBATCH --gpu-bind=single:1
#SBATCH --nodes=1
#SBATCH --exclusive

# Setup the job environment (this module needs to be loaded before any other modules)

module load PrgEnv-cray
module load rocm
module load craype-accel-amd-gfx90a
module load craype-x86-milan
module load cray-libsci

cd $SLURM_SUBMIT_DIR

export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=1
export OMP_PLACES=cores

# enable unified memory management
#export CRAY_ACC_USE_UNIFIED_MEM=1
#export HSA_XNACK=1

# ! TODO: add WORK_DIR
WORK_DIR=
PROJECT_DIR=$WORK_DIR/UCNS3D/PROTOTYPE/gpu_impl
EXECUTABLE=test
PATH_TO_EXECUTABLE=$PROJECT_DIR/$EXECUTABLE

# Input parameters
NUM_ELEMS=5000

srun --ntasks=4 --cpus-per-task=1 --hint=nomultithread --distribution=block:block $PATH_TO_EXECUTABLE $NUM_ELEMS
