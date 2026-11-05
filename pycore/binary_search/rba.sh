#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH --mem=1GB 
#SBATCH --time=10:00:00 
#SBATCH --account=cdm8
#SBATCH --partition=sla-prio
#SBATCH --output=%x_%j_slurm.out 
module load gams
python3 binary_search.py
