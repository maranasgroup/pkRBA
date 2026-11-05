#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH --mem=1GB 
#SBATCH --time=10:00:00 
#SBATCH --account=open
#SBATCH --partition=basic
#SBATCH --output=%x_%j_slurm.out 
python runRBA.py 
