# put all default model-specific parameters here (e.g., file paths). 
# Can be overridden later.

'''Naming conventions this script was designed for:
	- BiGG format generally followed unless otherwise specified
		- for stoichiometric models,
		- compartments (e.g., for metabolites/rxns) are after the last underscore in the ID. 
			- transport rxns include the compartment being transported from (when flux is positive) followed by the compartment being transported to (when flux is negative) (e.g., CO2t_c_e transports CO2 from cytosol to extracellular space)
'''

import os,sys
# add pycore to path
base_path = os.path.dirname(os.path.abspath(__file__))
if base_path not in sys.path:
	sys.path.append(base_path)
from cobrapy_custom_extras import load_cobra_model
from cobra import Model

class RBA_options:
	"""Stores default parameters for RBA model building and simulation.
	Attributes:
		biomass_impacted_by_PTMs (bool): Indicates if biomass is affected by post-translational modifications. Should be True for most models, but False is better for replicating behavior of RBA models before 2026.
		compartment_regex (dict): Dictionary of regular expressions for identifying compartments.
		exchange_regex (str): Regular expression pattern to identify exchange reactions.
		k_ribo_max (float): Upper limit during calculation of elongation rate for ribosomes, in amino acids per ribosome per hour. Helps calculate more realistic ribosome compositions and elongation rates. Doesn't need to be specified for most models, but could help if model-specific values are available. Default value derived from https://doi.org/10.1128/ecosal.5.2.3.
		organism_type (str): Type of organism (e.g., 'yeast', 'bacteria', 'mammal') to help identify ribosome subunits. Should be set to the type closest to your organism of interest, if none match.
		RBA_biomass_id (str): Identifier for the biomass reaction in the RBA model.
		SM_biomass_id (str): Identifier for the biomass reaction in the stoichiometric model.
		SM_path (str): Path to the stoichiometric model file.
		spont_GPRs (list): List of GPRs considered spontaneous.
		vmax (float): Default maximum reaction rate.
		remove_blocked_reactions (bool): Indicates whether to remove blocked reactions from the RBA model. Set to True or False.
		must_include_rxns (set): Set of reactions that must be included in the RBA model, even if they are blocked. Useful for including reactions that are important for the model but may not be active under certain conditions.
		always_localize_in_encoded_compartment (bool): Whether to always let proteins appear in the compartment where they're encoded. True is more lenient and less prone to arbitrarily excluding locations, but may be less accurate if only a protein's mRNA appears in its encoded compartment. False by default, since this can always be reversed during model building if desired.
		"""
	def __init__(self, SM_path, exchange_regex='^EX_', SM_biomass_id='BIOMASS', RBA_biomass_id='BIOMASS',
				 vmax=1000, compartment_regex={'default':'_[^_]+$'}, spont_GPRs=['SPONT', 'UNKNOWN', 'TRUE'], extra_proteins=[], protein_id_in_gene_format=None,
				 biomass_impacted_by_PTMs=True, k_ribo_max=22*3600, organism_type='', ribo_dict={}, cytosol='c', mito_compartments=[], decomp_compartments={}, remove_blocked_reactions=False, must_include_rxns=[], always_localize_in_encoded_compartment=False):
		# compartment_regex is a dictionary with keys for object types (e.g., 'metabolite', 'reaction', 'enzyme')
		## for any object type not specified, the default regex will be used
		# sets all input arguments as attributes
		for k, v in locals().items():
			if k != 'self':
				setattr(self, k, v)
	def load_SM(self):
		"""Load the stoichiometric model from the specified path."""
		model = load_cobra_model(self.SM_path)
		if not isinstance(model, Model):
			raise TypeError(f"Expected a cobra.Model instance, got {type(model)}")
		return model
	# check if model from SM_path has >2 compartments (e.g., more than cytosol and extracellular space)
	def has_subcellular_compartments(self):
		"""Check if the model has multiple compartments."""
		# self.protein_compartment = model.compartments[0]
		return len(model.compartments) > 2

# read '../build_model/input/options.yml' as dict to get the SM_path
options_path = os.path.join(base_path, '../build_model/input/options.yml')
import yaml
with open(options_path, 'r') as file:
	options = yaml.safe_load(file)
	SM_path = os.path.join(base_path, options.get('SM_path', None))  # Use default if not specified in options.yml
# check if SM_path exists
if not os.path.exists(SM_path):
	raise FileNotFoundError(f"Specified SM_path does not exist: {SM_path}")
# make an instance of the RBA_options class named model_defaults
RBA_defaults = RBA_options(SM_path=SM_path)
# assign any other options in options file to RBA_options attributes
for key, value in options.items():
	if key != 'SM_path':
		setattr(RBA_defaults, key, value)

main_script_globals = sys.modules['__main__'].__dict__
# check if 'load_model' key is set to False in the main script, to save time loading the model if not needed
load_model = main_script_globals.get('load_model', True)
if load_model:
	print(RBA_defaults.SM_path)
	model = RBA_defaults.load_SM()