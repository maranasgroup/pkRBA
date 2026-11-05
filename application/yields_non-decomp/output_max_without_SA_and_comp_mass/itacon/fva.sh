#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH --mem=1GB 
#SBATCH --time=60:00:00 
#SBATCH --account=open
#SBATCH --partition=basic
#SBATCH --output=%x_%j_slurm.out 
module load gams
gams FVA_max_yield.gms o=/dev/null --solver=soplex

