#!/bin/bash

#SBATCH -o Results_R1/Replicate_1_%a.Rout
#SBATCH --array=1-100
#SBATCH --mail-user=yiqiaozhang@ufl.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --account=k.ekvall
#SBATCH --mem-per-cpu=6gb
#SBATCH -t 48:00:00

module load R

R CMD BATCH --vanilla TimingComparison1_R1.R  Results_R1/Replicate_1_${SLURM_ARRAY_TASK_ID}.Rout
