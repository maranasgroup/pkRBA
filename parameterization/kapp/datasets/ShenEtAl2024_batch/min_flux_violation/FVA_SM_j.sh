#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH --mem=0
#SBATCH --time=10:00:00 
#SBATCH --account=cdm8
#SBATCH --partition=sla-prio
#SBATCH --output=%x_%j_slurm.out
gams FVA_SM_j.gms --solver=cplex

