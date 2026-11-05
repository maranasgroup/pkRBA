# Compare yields across various conditions
# add rxns making products of interest
import pandas as pd
import cobra
from copy import deepcopy
import json
from cobra.io import load_model
from pathlib import Path

import sys,os,shutil
sys.path.append('../')
# from pycore import *

# (optional) ensures output is visible when running as slurm job
from functools import partial
print = partial(print, flush=True)

# generate dataframe of yields from output_max_default, output_max_default_non-decomp, non-decomp_output_max_with_SA_no_comp_mass, non-decomp_output_max_without_SA_and_comp_mass
dirs = ('non-decomp_output_max_without_SA_and_comp_mass', 'non-decomp_output_max_with_SA_no_comp_mass_fMSA_min', 'output_max_default_non-decomp', 'output_max_default','non-decomp_output_max_without_SA_and_comp_mass_GSM_prot_levels')
overall = pd.DataFrame()
# arrange these so all their values neighbor each other
cols_to_combine=('yieldFBA','vCarbonSourcesRBA','rateRBA','yieldRBA')
datasets = {}
for d in dirs:
	datasets[d] = pd.read_csv(f"./{d}/compiled_results.csv")
	for c in datasets[d].columns:
		if c in cols_to_combine:
			datasets[d][f'{c}_{d}'] = datasets[d][c]
		if c not in overall.columns:
			overall[c] = datasets[d][c]
# sort so cols_to_combine are next to each other, in order of dirs
for c in cols_to_combine:
	for d in dirs:
		overall[f'{c}_{d}'] = datasets[d][f'{c}_{d}']
overall.set_index('product',inplace=True)				
overall.to_csv('yield_conditions_compared.tsv',sep='\t')