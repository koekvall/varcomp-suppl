#!/bin/bash

#SBATCH -o Run_R1_Zoomed_%a.Rout
#SBATCH --array=1-2
#SBATCH --mail-user=amolstad@ufl.edu
#SBATCH --mail-type=END
#SBATCH --account=amolstad
#SBATCH --qos=amolstad
#SBATCH --mem-per-cpu=2gb
#SBATCH -t 72:00:00

module load R

R CMD BATCH --vanilla CIsExample_R1_Zoomed.R  Run_R1_Zoomed_${SLURM_ARRAY_TASK_ID}.Rout
