#!/bin/bash -l
#SBATCH --job-name=multi-GPU-UCNS3D  # Job name
#SBATCH --output=multi-GPU-UCNS3D.o%j # Name of stdout output file
#SBATCH --error=multi-GPU-UCNS3D.e%j  # Name of stderr error file
#SBATCH --partition=standard-g
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --gpus-per-node=8
#SBATCH --time=01:00:00
#SBATCH --account=project_x

module load rocm

export LD_LIBRARY_PATH=/scratch/project_x/lib/:$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
export MPICH_GPU_SUPPORT_ENABLED=1
#export MPICH_ENV_DISPLAY=1
export OMP_TARGET_OFFLOAD=MANDATORY
export CRAY_ACC_USE_UNIFIED_MEM=0
export HSA_XNACK=0
export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_DYNAMIC=FALSE

cat > select_gpu <<'EOF'
#!/bin/bash
export ROCR_VISIBLE_DEVICES=$SLURM_LOCALID
export OMP_DEFAULT_DEVICE=0
exec "$@"
EOF
chmod +x ./select_gpu

CPU_BIND="map_cpu:49,57,17,25,1,9,33,41"


##--gres-flags=allow-task-sharing
srun --cpu-bind=${CPU_BIND} \
     ./select_gpu ./ucns3d_p
