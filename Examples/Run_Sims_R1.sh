#!/bin/bash

#SBATCH -o Run_R1_%a.Rout
#SBATCH --array=1-2
#SBATCH --mail-user=amolstad@ufl.edu
#SBATCH --mail-type=END
#SBATCH --account=amolstad
#SBATCH --qos=amolstad
#SBATCH --mem-per-cpu=2gb
#SBATCH -t 72:00:00

module load R

R CMD BATCH --vanilla CIsExample_R1.R  Run_R1_${SLURM_ARRAY_TASK_ID}.Rout
