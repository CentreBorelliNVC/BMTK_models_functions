#!/bin/bash
#SBATCH --job-name=run_pointnet_job
#SBATCH --account=MouseCortexModel
#SBATCH --mail-type=ALL
#SBATCH --mail-user=margaux.calice@etu.u-paris.fr
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --partition=cpu_short
#SBATCH --output=output_run_pointnet.txt
#SBATCH --error=output_error_run_pointnet.txt

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

cd ${SLURM_SUBMIT_DIR}

conda activate ENV_NEST33
srun run_pointnet.py
