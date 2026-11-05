"""Custom functions extending COBRApy's functionality."""

from cobra import Model

def save_cobra_model(model:Model, path:str=None, *args, **kwargs) -> None:
	"""Save a stoichiometric model using COBRApy. If no path given, defaults to the model's ID followed by .json."""
	import cobra
	# check if input is a string
	if path is None:
		path = str(model.id) + '.json'
	if isinstance(path, str):
		# check file extension to determine file format
		ext = path.split('.')[-1] if '.' in path else 'json'
		extension_mapping = {'json': cobra.io.save_json_model, 
					   'xml': cobra.io.write_sbml_model, 
					   'mat': cobra.io.save_matlab_model,
					   'yaml': cobra.io.save_yaml_model}
		if ext in extension_mapping:
			extension_mapping[ext](model, path, *args, **kwargs)
		else:
			raise ValueError('File format not recognized')
	else:
		raise ValueError('Input must be a file path string')

def load_cobra_model(self, *args, **kwargs) -> Model:
	"""Load a stoichiometric model using COBRApy, either as a dict or a string with a file path."""
	import cobra
	# check if input is a string
	if isinstance(self, str):
		# check file extension to determine file format
		ext = self.split('.')[-1]
		extension_mapping = {'json': cobra.io.load_json_model, 
					   'xml': cobra.io.read_sbml_model, 
					   'sbml': cobra.io.read_sbml_model, 
					   'mat': cobra.io.load_matlab_model, 
					   'yaml': cobra.io.load_yaml_model}
		if ext in extension_mapping:
			model = extension_mapping[ext](self, *args, **kwargs)
		else:
			raise ValueError('File format not recognized')
	elif isinstance(self, dict):
		model = cobra.io.dict_to_model(self)
	else:
		raise ValueError('Input must be a dictionary or a file path string')
	return model

def rename_gene(model, old_id, new_id=None):
	"""Rename a gene in a COBRApy model (in place), with optional new ID. If new_id not given, will attempt to find a shorter name from gene attributes."""
	from cobra.manipulation import rename_genes
	if old_id not in model.genes:
		print(f"Warning: gene {old_id} not found in model; cannot rename to {new_id}")
		return
	gene = model.genes.get_by_id(old_id)
	if new_id is None:
		for attr in [gene.name, getattr(gene, 'annotation', {}).get('ncbiprotein', '')]:
			if attr and len(attr) < len(old_id) and attr not in model.genes:
				new_id = attr
				break
		new_id = gene.name if len(gene.name) < len(old_id) else old_id
	if new_id in model.genes:
		print(f"Warning: gene {new_id} already exists in model; cannot rename {old_id} to {new_id}")
		return
	# move old_id to notes for reference, in case of future issues with gene ID changes
	if not hasattr(gene,"notes"):
		gene.notes = {}
	gene.notes.setdefault("notes","")
	msg = f"Original gene ID was {old_id}"
	if msg not in gene.notes['notes']:
		gene.notes['notes'] += msg
	rename_genes(model, {old_id: new_id})