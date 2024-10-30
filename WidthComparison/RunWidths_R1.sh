#!/bin/bash

#SBATCH -o Results_R1/Width_%a.Rout
#SBATCH --array=1-40
#SBATCH --mail-user=amolstad@ufl.edu
#SBATCH --mail-type=END
#SBATCH --account=amolstad
#SBATCH --qos=amolstad-b
#SBATCH --mem-per-cpu=5gb
#SBATCH -t 96:00:00

module load R

R CMD BATCH --vanilla Widths_R1.R  Results_R1/Width_${SLURM_ARRAY_TASK_ID}.Rout
