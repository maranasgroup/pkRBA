#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH --mem=0
#SBATCH --account=cdm8
#SBATCH --partition=sla-prio
#SBATCH --time=08:00:00 
#SBATCH --output=%x_%j_slurm.out
python3 B2_min_flux_violation.py && python3 C1_calculate_kapp.py
