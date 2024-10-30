#!/bin/bash

#SBATCH -o  Results_R1/Replicate_%a.Rout
#SBATCH --array=1-144
#SBATCH --mail-user=amolstad@ufl.edu
#SBATCH --mail-type=END
#SBATCH --account=amolstad
#SBATCH --qos=amolstad-b
#SBATCH --mem-per-cpu=5gb
#SBATCH -t 96:00:00

module load R

R CMD BATCH --vanilla Eigen_Sim1_R1.R  Results_R1/Replicate_${SLURM_ARRAY_TASK_ID}.Rout
