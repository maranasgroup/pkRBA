# Model-specific settings for kapp calculations, as well as steps that must be run each time
import pandas as pd
import numpy as np
import requests,json,os,shutil,sys

# retrieve name of the directory containing the options file
options_dir = os.path.dirname(os.path.abspath(__file__)) + '/'
model_root_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../..')) + '/' # add trailing slash for easier path concatenation
path_model = model_root_path+'GAMS/model/'
pycore_path = model_root_path+'pycore/'
for path in [model_root_path, pycore_path]:
    if path not in sys.path:
        sys.path.append(path)

from pycore import * 

vmax = RBA_defaults.vmax # max flux in either direction, before applying nscale
nscale=1e4 # will be used as scale when running GAMS
nscale2=nscale # 2nd scale, if different one needed; otherwise, set to same as nscale

# Load path
path_gen = model_root_path+'build_model/'
path_gams = model_root_path+'GAMS/'

prot_path = path_gen + 'model/PROTEIN_stoich_curation.tsv'
unknown_prot_path = path_gen + 'input/PROTEIN_stoich_curation_unknown.xlsx'
model_xlsx_path = path_model + 'RBA_stoichiometry.tsv'
metab_rxns_path = path_model + 'RBA_rxns_rxnmetabolicnetwork.txt'
gsm_rxn_ids_path = path_model + 'SM_rxn_ids.txt'
sij_path = path_model + 'RBA_sij.txt'
aa_mapping_path = path_gen + 'input/PROTEIN_amino_acid_map.txt'
aa_map = pd.read_csv(aa_mapping_path, sep='\t')
aa_dict = dict(zip(aa_map['aa_abbv'], aa_map['MW']))

# for flux data; keep path blank if none provided
flux_path = os.path.abspath(os.path.join(options_dir, '../raw_data_files/fluxes_net.tsv'))
flux_index = 'rxn_id'
## upper and lower estimates on measured flux. Set to the same column if only one set of flux estimates provided, or to different columns if using multiple estimates (e.g., if confidence intervals or a range of values are provided).
col_v_exp_lo = 'val'
col_v_exp_up = 'val'
## upper and lower bounds on fluxes (e.g., from model used in MFA), to distinguish actual fluxes from errors (e.g., in solver precision or measurement)
col_LB = 'vLB'
col_UB = 'vUB'

# must match the growth rate in other files 
# 	may help to add or subtract measurement error in the value itself, for future reference (e.g., mu = 0.060 - 0.001, not mu = 0.059)
# if/else statements below help automatically determine mu; they set mu to 1st col_v_exp_lo value of flux_path where flux_index == 'BIOMASS', in case multiple growth rates are provided
if flux_path != '' and os.path.exists(flux_path):
    df_flux = pd.read_csv(flux_path, sep='\t')
    mu = df_flux.loc[df_flux[flux_index]=='BIOMASS', col_v_exp_lo].values[0]
else: mu = 0.1
print(f"Using mu = {mu} 1/hr based on flux data file {flux_path}")
# protein fraction (disable by uncommenting "ptot = 1" unless composition varies w/ growth rate)
ptot = 0.484 # from MFA biomass
# ptot = 0.524 # from measured biomass

# protein data
df_raw = pd.read_excel(os.path.abspath(os.path.join(options_dir, '../raw_data_files/ShenEtAl2024_multiomics_41589_2024_1571_MOESM3_ESM.xlsx')),sheet_name='Table 10c. abs_prot_IO_SD108')
# make sure df_raw.index maps to either the 'id' column or the 'uniprot' column at prot_path
df_raw.index = df_raw['Entry'].to_list() # name of protein
# df_raw.columns = df_raw.iloc[0]
cols_data = ['replicate 1'] # where protein abundance data is stored
# cols_data = ['mean (renormalized)'] # where protein abundance data is stored
uniprot_col = 'Entry' # set to '' if no column with accession names provided
proteomics_units = 'g/g protein' # supported options: 'g/g protein' (originally built with this in mind), 'g/gDW'
use_ribo_data = False # False to ignore ribosome subunit abundance data (still used for MPFS constraints)
find_kribo = True # if true, runs kribo bisection

# if True, proteins with measured translation rates of 0 can be made at a flux equal to prosynSlackUB, thus discouraging their use where possible. If you want such proteins to be made, you thus may have to specify them in the phenotype.txt file as having fixed values for prosynSlackUB. 
# if False, proteins with measured translation rates of 0 can't be made. This may cause issues if no isozymes exist for such proteins but their rxns are essential, and such proteins may still be present but unmeasured due to measurement errors.
allow_trans_when_measurement_is_0 = True 

nonmodeled_proteome_allocation = 0
# nonmodeled_proteome_allocation = 0.503973558
dummy_protein = {'id':'PROSYN-PROTDUMMY','AA abundances':dict()}

path_data = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), './protein_data.tsv'))
if os.path.exists(path_data):
    df_data = read_spreadsheet(path_data)
    df_data.index = df_data['id'].to_list()
    df_data = df_data[df_data['conc (g/gDW)'] > 0]
    # Excluding ribosome protein subunit (conflicting if fit to both enzymatic and ribosomal protein data)
    # make perfect copy of df_data
    df_data_full = df_data.copy()
    if not use_ribo_data:
        df_data = df_data[(df_data.type == 'truedata_enz') | (df_data.type == 'gapfill_subunit')]

biom_id = RBA_defaults.RBA_biomass_id

# set to True to avoid printing .lst files and other large log files
hide_output = True
#doesn't work for GAMS 42 (on Roar cluster): output_redirect_str = ' writeOutput=0 ll=0' if hide_output else ''
#output_redirect_str = ' o=/dev/null ll=0 lo=0' if hide_output else ''
output_redirect_str = ' o=/dev/null' if hide_output else ''

spont_rxn_suffixes = RBA_defaults.spont_GPRs
rxns_to_ignore_for_kapps = [] # ignored when printing rxns during C1
