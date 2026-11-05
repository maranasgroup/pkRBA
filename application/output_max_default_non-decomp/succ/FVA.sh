#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH --mem=1GB 
#SBATCH --time=30:00:00 
#SBATCH --account=open
#SBATCH --partition=basic
#SBATCH --output=%x_%j_slurm.out 
module load gams
gams FVA_max_yield.gms --mu=0.143749481
