# add rxns making products of interest
import pandas as pd
import cobra
from copy import deepcopy
import json
from cobra.io import load_model
from pathlib import Path

import sys,os,shutil
model_root_path = Path(os.path.abspath('../../../'))
sys.path.append(model_root_path.as_posix())
from pycore import *

from RBA_defaults_from_FBA import *

# use find_pathways function to make a graph of the pathways making the product of interest (the rxn in question is vprod_coreid), and find the rxns that are in the pathways
modelname = os.path.splitext(os.path.basename(SM_path))[0]
newmodel_path = '../../input/'+modelname+'_pathwayadd.json'
model2 = load_cobra_model(newmodel_path)

# load fluxes from RBA result
with open('runRBA_max_prod.flux.txt') as f:
	text = f.read().split('\n')
	# 1st column is rxn id, 3rd column is flux value
	text = [i.split('\t') for i in text if i != '']
net_flux = {}
gross_flux = {}
for i in text:
	_, base_id, dir, _ = extract_details_from_rxnid(i[0])
	if base_id not in model2.reactions:
		continue
	net_flux.setdefault(base_id,0)
	net_flux[base_id] += float(i[2]) if dir == 'FWD' else -float(i[2])
	gross_flux.setdefault(base_id,{'FWD':0,'REV':0})
	gross_flux[base_id][dir] += float(i[2])

prod_in_model = [m.id for m in model2.reactions.get_by_id(vprod_coreid).metabolites][0]
target_reactants = []
for c in c_sources_list:
	for m in model2.reactions.get_by_id(c['rxn']).metabolites:
		target_reactants.append(m.id)

print('Finding pathways to produce',prod_in_model)

node_font_size=8
edge_font_size=8
exclude_common_mets=True
# load /Users/ejm6426/Documents/pkRBA/application/non-decomp_output_max_with_SA_no_comp_mass/3hppa/pathways-making-3hpp_e-in-iIsor850-G3.gml, exclude all edges with identical labels and fluxes (rounded to 6 decimal places), and save to /Users/ejm6426/Documents/pkRBA/application/non-decomp_output_max_with_SA_no_comp_mass/3hppa/pathways-making-3hpp_e-in-iIsor850-G3_excluded-identical-fluxes.gml
import networkx as nx

mg = find_pathways(model2,[prod_in_model],flux_dict=net_flux,target_reactants=target_reactants,flux_reporting_decimals=4,exclude_common_mets=exclude_common_mets,
				   node_font_size=node_font_size,edge_font_size=edge_font_size,output_path='./')

# color decompartmentalized rxns in the graph
key_rxns_path = '../../input/GAMS_model_application/decompartmentalized_rxns.txt'
if os.path.exists(key_rxns_path):
	G_chosen = mg[-1]
	with open (key_rxns_path) as f:
		text = [i.strip("'") for i in f.read().split('\n') if i not in ('','/\n','/')]
	key_rxns = set()
	key_color = '#FF0000'
	for i in text:
		tag, base_id, dir, _ = extract_details_from_rxnid(i)
		if tag in ('RXN', 'RXNADD') and base_id in model2.reactions:
			if gross_flux.get(base_id, {'FWD': 0, 'REV': 0})[dir] != 0:
				key_rxns.add(base_id)

	for u, v, key, data in G_chosen.edges(keys=True, data=True):
		rxn_id = data.get('id', '')
		if rxn_id in key_rxns:
			data['graphics'] = {'fill': key_color, 'width': 2.0, 'targetArrow': 'standard', 'sourceArrow': 'none'}
			# set label color to key_color
			if 'label' in data:
				data.setdefault('LabelGraphics', {}).setdefault('color', key_color)

	# output the graph to a GML file
	output_gml_path = f'pathways-making-{prod_in_model}-in-{model2.id}.gml'
	nx.write_gml(G_chosen, output_gml_path)
# Optional: specify a path to an existing GML file to load and compare fluxes
gml_path = None
# example path: 
gml_path = f'{model_root_path}/application/non-decomp_output_max_with_SA_no_comp_mass/3hppa/pathways-making-3hpp_e-in-iIsor850-G3.gml'
if gml_path is not None:
	G = nx.read_gml(gml_path)
	alt_flux_dict = {}

	# Create a new graph with the same nodes
	G_new = nx.DiGraph()
	G_new.add_nodes_from(G.nodes(data=True))

	# Add edges, excluding those with identical labels and fluxes (rounded to n decimal places)
	decimal_places = 1
	for u, v, data in G.edges(data=True):
		rxn_id = data.get('id', '')
		flux = data.get('flux', 0)
		alt_flux_dict[rxn_id] = flux
		if model2.reactions.has_id(rxn_id) and round(flux, decimal_places) == round(net_flux.get(rxn_id, 0), decimal_places):
			continue
		G_new.add_edge(u, v, **data)
	# remove any nodes that have no edges
	nodes_to_remove = [node for node in G_new.nodes() if G_new.degree(node) == 0]
	G_new.remove_nodes_from(nodes_to_remove)

	# Save the new graph to a GML file
	nx.write_gml(G_new, 'pathways-making-3hpp_e-in-iIsor850-G3_excluded-identical-fluxes.gml')
	# make a version of this graph that also excludes edges with identical labels and fluxes (rounded to n decimal places), and save to this folder
	G_new = nx.DiGraph()
	pathway_graphs = mg if isinstance(mg, list) else [mg]
	for pathway_graph in pathway_graphs:
		if not hasattr(pathway_graph, 'nodes'):
			continue
		G_new.add_nodes_from(pathway_graph.nodes(data=True))
		for u, v, data in pathway_graph.edges(data=True):
			rxn_id = data.get('id', '')
			flux = data.get('flux', 0)
			if round(flux, decimal_places) == round(alt_flux_dict.get(rxn_id, 0), decimal_places):
				continue
			print(f'Adding edge {u} -> {v} with flux {flux} (rounded: {round(flux, decimal_places)})')
			G_new.add_edge(u, v, **data)
	# remove any nodes that have no edges
	nodes_to_remove = [node for node in G_new.nodes() if G_new.degree(node) == 0]
	G_new.remove_nodes_from(nodes_to_remove)
	nx.write_gml(G_new, 'pathways-making-3hpp_e-in-iIsor850-G3_excluded-identical-fluxes2.gml')