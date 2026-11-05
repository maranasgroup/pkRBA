# updates model-specific settings in kapp_options.py, then runs kribo minimization (if user chooses to)
load_model = False
from kapp_options import *

path_rxns_list = path_gams + 'model/RBA_rxns.txt'
rxns_with_enz = [i.strip("'") for i in open(path_model+'RBA_rxns_with_enz.txt').read().split('\n') if i not in ('','/')]
# read protein length as dict for use in kribo calculation, separating key from value using space, and removing quotes from keys
protein_length = dict()
with open(path_gams + 'model/RBA_proteinLength.txt') as f:
    for i in f.read().split('\n'):
        if i not in ['','/']:
            k,v = i.split(' ')
            protein_length[k[1:-1].split('-',1)[1]] = float(v)
pro_trans_dict = dict() # for assigning proteins to ribosome types based on translation reactions, for use in kribo calculation
with open(path_gams + 'model/translation.txt') as f:
    for i in f.read().split('\n'):
        if i not in ['','/']:
            k,v = i.split("'.'")
            pro_trans_dict.setdefault(v[:-1].split('-',1)[1], []).append(k[1:])
with open(sij_path, 'r') as f:
    # make sij_dict, a dict where each key is the rxn id (after "'.'") and each value is a dict with keys 'met_id' and 'sij' for the corresponding met_id and sij value for that rxn-met pair, to be used in kribo calculation
    sij_dict = dict()
    for line in f.read().split('\n'):
        if line not in ['','/']:
            met_id, rxn_id, sij = line.replace("'.'", "' '").split(' ')
            sij_dict.setdefault(rxn_id[1:-1], {})[met_id[1:-1]] = float(sij)

path_out = './min_flux_violation/'
path_enz_level = './enz_from_proteome/enz_flux_calculation.txt'

run_setting_file_from = './GAMS_setting_files/min_flux_violation_GAMS_settings.txt'
run_setting_file_to = './min_flux_violation/min_flux_violation_GAMS_settings.txt'

df_prot = read_spreadsheet(prot_path)
df_prot.index = df_prot.id.to_list()
min_kapps = dict() # for determining allocation of non-measured proteins w/ used rxns

#### Create directory and copy run settings
if os.path.isdir(path_out) == False:
    os.makedirs(path_out)
shutil.copy(run_setting_file_from, run_setting_file_to);

data_val = []; data_idx = [];
for i in df_data.index:
    data_idx.append("'PROSYN-" + df_data.id[i] + "'")
    data_val.append("'PROSYN-" + df_data.id[i] + "' " + str(df_data.loc[i, 'vtrans (mmol/gDW/h)']))

data_val = ['/'] + data_val + ['/']
data_idx = ['/'] + data_idx + ['/']
with open(os.path.join(path_out, 'proteome_data.txt'), 'w') as f:
    f.write('\n'.join(data_val))
with open(os.path.join(path_out, 'proteome_data_set.txt'), 'w') as f:
    f.write('\n'.join(data_idx))

#### Determine likely active and inactive reactions informed by proteomics data and calculation in B1_enz_from_proteome.py
# List out model reactions
with open(path_rxns_list) as f:
    idx = f.read().split('\n')[1:-1]
idx = [i[1:-1] for i in idx]
rxns_all = [i for i in idx if i.split('-')[0] == 'RXN']
rxns_enz_all = [i for i in idx if i.split('-')[0] == 'ENZLOAD']
enz_rxn_pairs = dict()
# ENZLOAD reactions are enzyme-specific, whereas RXN reactions are not.  Keep
# an explicit mapping rather than inferring an enzyme from an RXN identifier.
rxn_to_enzloads = dict()
enzload_to_rxn = dict()
# define all enz-rxn pairs, for determining prototype kapps
for enzload in rxns_enz_all:
    enzload_id = enzload.split('-', 1)[1]
    rxn_id, enz_id = enzload_id.rsplit('-', 1)
    rxn = 'RXN-' + rxn_id
    tag,rxn_base_id,rxn_dir,_ = extract_details_from_rxnid(enzload)
    rxn_to_enzloads.setdefault(rxn, []).append(enzload)
    enzload_to_rxn[enzload] = rxn
    # find all rxns catalyzed by the same enzyme
    if enz_id not in enz_rxn_pairs.keys():
        enz_rxn_pairs[enz_id] = dict()
    enz_rxn_pairs[enz_id][enzload] = {'rxn':rxn,'tag':tag,'rxn_base_id':rxn_base_id,'rxn_dir':rxn_dir,'v':0,'v_enzload':0,'kapp_proto':0}

rxns_enz_rxns = {pair['rxn'] for pairs in enz_rxn_pairs.values() for pair in pairs.values()}
rxns_nonenz = [i for i in rxns_all if i not in rxns_enz_rxns]

# List active enzymatic reactions (rxns with enzymes that can be made, using estimates from B1_enz_from_proteome.py)
with open(path_enz_level) as f:
    enz_fluxes = f.read().split('\n')
enz_fluxes = [i for i in enz_fluxes if i != '']

rxns_enz_active = dict()
for i in enz_fluxes:
    enzid,v = i.split('\t')
    etype = enzid.split('-')[0]
    if etype == 'ENZLOAD':
        rxn = enzload_to_rxn[enzid]
        enz_id = enzid.split('-', 1)[1].rsplit('-', 1)[1]
        if float(v) > 0:
            rxns_enz_active[enzid] = {'enz':enz_id,'rxn':rxn.rsplit('-',1)[0],'enzload':enzid,'v_enzload':float(v),'v':0}
            enz_rxn_pairs[enz_id][enzid]['v_enzload'] = float(v)

# List of enzymatic reactions with no proteomics data 
rxns_enz_with_no_prodata = [i for i in rxns_all if i not in [x['rxn'] for x in rxns_enz_active.values()] + rxns_nonenz]
enzload_with_no_prodata = ["'" + e + "'" for r in rxns_enz_with_no_prodata for e in rxn_to_enzloads.get(r, [])]
rxns_enz_with_no_prodata = ["'" + i + "'" for i in rxns_enz_with_no_prodata]
rxns_enz_with_no_prodata = ['/'] + rxns_enz_with_no_prodata + ['/']
fname = os.path.join(path_out, 'enzload_enz_with_no_prodata.txt')
with open(fname, 'w') as f:
    f.write('\n'.join(['/'] + enzload_with_no_prodata + ['/']))
fname = os.path.join(path_out, 'rxns_enz_with_no_prodata.txt')
with open(fname, 'w') as f:
    f.write('\n'.join(rxns_enz_with_no_prodata))

fname = './min_flux_violation/min_flux_violation.modelStat.txt'
# remove to prevent past runs from overriding status checks
if os.path.isfile(fname):
    os.remove(fname)

cmds = ['cd ' + path_out,
        'module load gams',
        'gams min_flux_violation.gms --nscale=' + str(nscale) + output_redirect_str]
os.system('\n'.join(cmds))
stop_if_run_failed('./min_flux_violation/min_flux_violation.modelStat.txt')
#### Convert GAMS-scaled flux to actual flux
# All fluxes
rxns_used = []
rxns_enz_with_no_prodata_essential = dict()
enzyme_rxns = []
enz_rxn_dict = {}

fname = os.path.join(path_out, 'min_flux_violation.rxns_essential_with_no_prodata_gamsscaled.txt')
f2name = os.path.join(path_out, 'min_flux_violation.enzload_essential_with_no_prodata_gamsscaled.txt')
with open(fname) as f, open(f2name,'w') as f2:
    newlines = []
    for i in f.read().split('\n'):
        if i not in ['','/']:
            j,v = i.split('\t')
            rxns_enz_with_no_prodata_essential[j[1:-1]] = float(v) 
            for enzload in rxn_to_enzloads.get(j[1:-1], []):
                newlines.append("'" + enzload + "'\t" + v)
        else:
            newlines.append(i)
    f2.write('\n'.join(newlines))
    # print(rxns_enz_with_no_prodata_essential)

fname = os.path.join(path_out, 'min_flux_violation.flux_gamsscaled.txt')
with open(fname) as f:
    fluxes = f.read().split('\n')
fluxes = [i for i in fluxes if i != '']
fluxes_new = []
for i in fluxes:
    r,vtype,val = i.split('\t')
    val = float(val) / nscale
    fluxes_new.append('\t'.join([r, vtype, str(val)]))
    rxns_used.append("'" + r + "'")
    # if rxn (as determined by characters after last "-") isn't SPONT or UNKNOWN
    if r in rxns_with_enz:
        for enz,enz_details in enz_rxn_pairs.items():
            for enzload,details in enz_details.items(): 
                if details['rxn'] != r: continue
                enz_rxn_dict.setdefault(enz, []).append(enzload)
                enz_rxn_pairs[enz][enzload]['v'] = float(val)
                if enzload in rxns_enz_active:
                    rxns_enz_active[enzload]['v'] = val
                enzyme_rxns.append("'" + enzload + "'")
# find min kapps to construct constraints that get the essential enz rxns w/ no prodata to have kapps as close as possible to them
for j in rxns_enz_active.keys():
    v = rxns_enz_active[j]['v']
    if v > 0:
        min_kapp = mu * v / (float(rxns_enz_active[j]['v_enzload']))
        enz_rxn_pairs[rxns_enz_active[j]['enz']][j]['kapp_proto'] = min_kapp
        min_kapps[j] = min_kapp
kapp_benchmark = float(np.median(list(min_kapps.values())))
kapp_max = float(max(list(min_kapps.values())))
print('max predicted kapp',str(kapp_max))
print('benchmark kapp',str(kapp_benchmark))

enzload_constraints = [] # determining enzyme allocation for next problem
kapp_init = []
c = 0

from copy import deepcopy
enz_rxn_pairs_backup = deepcopy(enz_rxn_pairs)
# find kapps for all rxns that don't have any
for enz in enz_rxn_pairs.keys():
    for rxn in enz_rxn_pairs[enz].keys():
        if enz_rxn_pairs[enz][rxn]['kapp_proto'] == 0:
            rxn_base_id = enz_rxn_pairs[enz][rxn]['rxn_base_id']
            v = enz_rxn_pairs[enz][rxn]['v']
            v_enzload = enz_rxn_pairs[enz][rxn]['v_enzload']
            kapps = [enz_rxn_pairs_backup[enz][i]['kapp_proto'] for i in enz_rxn_pairs[enz].keys() if enz_rxn_pairs_backup[enz][i]['kapp_proto'] > 0]
            # assign kapp as the sum of all fluxes catalyzed by the enzyme, divided by the enzyme's synthesis flux
            for enz2 in enz_rxn_pairs.keys():
                for rxn2 in enz_rxn_pairs[enz2].keys():
                    if enz_rxn_pairs[enz2][rxn2]['rxn_base_id'] == rxn_base_id:
                        if enz_rxn_pairs_backup[enz2][rxn2]['kapp_proto'] > 0:
                            kapps.append(enz_rxn_pairs_backup[enz2][rxn2]['kapp_proto'])
                            break
            # assign benchmark kapp if no other kapp values are available
            enz_rxn_pairs[enz][rxn]['kapp_proto'] = min([kapp_benchmark] + kapps)
        # constraint to encourage equal distribution of enzload amongst all used rxns
        c += 1
        enzsyns = [r for r in sij_dict.keys() if 'ENZ-'+enz in sij_dict[r].keys() and r.startswith('ENZSYN')]
        enzload = rxn
        enzload_constraints.append(f"Equation EnzLoadConstraint{c}; EnzLoadConstraint{c}.. EnzLoadSlackPos('{enzload}') - EnzLoadSlackNeg('{enzload}') + v('{enzload}') =e= (" + '+'.join([f"v('{enzsyn}')" for enzsyn in enzsyns]) + f") / {len(enz_rxn_pairs[enz].keys())};")
        kapp_init += [f"'{enzload}' {str(enz_rxn_pairs[enz][rxn]['kapp_proto'])}"]

# for enz,rxns in enz_rxn_dict.items():
#     enzsyn = 'ENZSYN-' + enz
#     for rxn in rxns:
#         c += 1
#         c1 += 1
#         eload = rxn.replace('RXN-','ENZLOAD-',1)
#         j = rxn.replace('ENZLOAD-','RXN-',1)
#         if j in rxns_enz_with_no_prodata_essential.keys():
#             kapp_init.append("Equation kappEst" + str(c1) + "; kappEst" + str(c1) + ".. kappEstSlackPos('" + j + "') - kappEstSlackNeg('" + j + "') + v('" + eload + "') =e= v('" + j + "') / (%mu% * " + str(kapp_benchmark) + ");") 
#         elif rxns_enz_active[j]['v'] > 0:
#             kapp_init.append("Equation kappEst" + str(c1) + "; kappEst" + str(c1) + ".. kappEstSlackPos('" + j + "') - kappEstSlackNeg('" + j + "') + v('" + eload + "') =e= v('" + j + "') / (%mu% * " + str(min_kapps[j]) + ");") 
#         # constraint to encourage equal distribution of enzload amongst all used rxns
#         enzload_constraints.append("Equation EnzLoadConstraint" + str(c) + "; EnzLoadConstraint" + str(c) + ".. " + "EnzLoadSlackPos('" + rxn + "') - EnzLoadSlackNeg('" + rxn + "') + v('" + rxn + "') =e= v('" + enzsyn + "') / " + str(len(rxns)) + ";")
fname = os.path.join(path_out, 'kapp_init.txt')
with open(fname, 'w') as f:
    f.write('\n'.join(['/'] + kapp_init + ['/']))
fname = os.path.join(path_out, 'enz_alloc_constraints.txt')
with open(fname, 'w') as f:
    f.write('\n'.join(enzload_constraints))
fname = os.path.join(path_out, 'min_flux_violation.flux.txt')
with open(fname, 'w') as f:
    f.write('\n'.join(fluxes_new))
fname = os.path.join(path_out, 'min_flux_violation.rxns_used.txt')
with open(fname, 'w') as f:
    f.write('\n'.join(['/'] + rxns_used + ['/']))
# for determining enzyme allocation
fname = os.path.join(path_out, 'min_flux_violation.enzyme_rxns.txt')
with open(fname, 'w') as f:
    f.write('\n'.join(['/'] + enzyme_rxns + ['/']))
   
# Only flux of reactions not associated with expressed enzymes
fname = os.path.join(path_out, 'min_flux_violation.flux_essential_with_no_prodata_gamsscaled.txt')
with open(fname) as f:
    fluxes = f.read().split('\n')
fluxes = [i for i in fluxes if i != '']
fluxes_new = []
for i in fluxes:
    r,vtype,val = i.split('\t')
    val = float(val) / nscale
    fluxes_new.append('\t'.join([r, vtype, str(val)]))
fname = os.path.join(path_out, 'min_flux_violation.flux_essential_with_no_prodata.txt')
with open(fname, 'w') as f:
    f.write('\n'.join(fluxes_new))
# for determining enzyme allocation
fname = './min_flux_violation/enz_alloc.modelStat.txt'
# remove to prevent past runs from overriding status checks
if os.path.isfile(fname):
    os.remove(fname)

cmds = ['cd ' + path_out,
        'module load gams',
        f'gams enz_alloc.gms --nscale={str(nscale)}{output_redirect_str}']
os.system('\n'.join(cmds))
stop_if_run_failed(fname)
# All fluxes
fname = os.path.join(path_out, 'enz_alloc.flux_gamsscaled.txt')
rxns_used = []
enzyme_rxns = []
enz_rxn_dict = {}
flux_dict = dict()
with open(fname) as f:
    fluxes = f.read().split('\n')
    flux_dict = {i.split('\t')[0]:float(i.split('\t')[2])/nscale2 for i in fluxes if i not in ['','/']}
V_ribo_K_ribo = {k:{'min': 0, 'max': 0} for k in set([i for sublist in pro_trans_dict.values() for i in sublist])} # ribosome flux * elongation rate for each ribosome type, to be used in kribo calculation
ribo_subunit_candidates = {k:{'SSU': {}, 'LSU': {}, 'unknown': {}, 'composition_dict': {}} for k in V_ribo_K_ribo.keys()} # list of proteins that could be subunits for each ribosome type, for use in kribo calculation
fluxes = [i for i in fluxes if i != '']
fluxes_new = []
aa_comp = dict()
biomass_from_RBA = dict()
ribo_name_conversions = dict() # for converting between ribosome naming schemes
ribo_name_conversions_preferred = dict() # for converting between ribosome naming schemes, with a preferred name for each set of conversions (for use in figures and interpretation)
ribo_name_organism_conversions = dict() # for converting ribosome names based on organism (e.g., mito vs cyto), for use in interpreting proteomics data and assigning ribosome subunit membership, where needed
organism_types = set()

def interpret_subunit_name(name, return_format='dict'):
    '''Split ribosome subunit names (using format from https://bangroup.ethz.ch/research/nomenclature-of-ribosomal-proteins.html) into their constituent parts (e.g., organism/organelle, large/small subunit, designated number, mito/chloroplast version, paralog indicated for plants)'''
    universal_ribo_regex = r'^(.)(.)(\d+)([cm]{0,1})(\w{0,1})$'
    origin, subunit_group, number, mito_chloro, paralog = re.match(universal_ribo_regex, name).groups() if re.match(universal_ribo_regex, name) else (None, None, None, None, None)
    subunit_dict = {'origin': origin, 'subunit_group': subunit_group, 'number': number, 'mito_chloro': mito_chloro, 'paralog': paralog}
    if return_format == 'dict':
        return subunit_dict
    if return_format == 'string, paralog-agnostic':
        return_str = ''.join(i if i not in [None, ''] else '' for i in [origin, subunit_group, number])
        return_str = return_str if return_str != '' else name
        return return_str

ribo_name_df = pd.read_csv(model_root_path + 'build_model/ribosome_name_conversions_non-plant.csv')
ribo_name_df = ribo_name_df.fillna('')
# convert into a dict where each cell (outside of "new name" column) is a key and its value is its corresponding "new name" value, for use in converting ribosome names in the proteomics data to the names used in the model
for i in ribo_name_df.index:
    new_name = re.sub(r'\s*\([^)]*\)', '', ribo_name_df.loc[i, 'new name'])
    for col in ribo_name_df.columns:
        if col != 'new name' and ribo_name_df.loc[i, col] != '':
            for old_name in re.split(r',\s*', ribo_name_df.loc[i, col]):
                if old_name in [new_name, '']:
                    continue
                old_name_clean = re.sub(r'\s*\([^)]*\)', '', old_name)
                organism_type = col.split()[0].lower() if len(col.split()) > 0 else None
                organism_types.add(organism_type)
                if organism_type == RBA_defaults.organism_type.lower():
                    if new_name in ribo_name_conversions_preferred.keys() and ribo_name_conversions_preferred[new_name] != new_name:
                        print(f"Warning: {new_name} is associated with multiple preferred names: {ribo_name_conversions_preferred[new_name]} and {new_name}. Keeping the last one.")
                    ribo_name_conversions_preferred[old_name_clean] = new_name
                ribo_name_conversions[old_name_clean] = new_name
                if ribo_name_organism_conversions.get(organism_type, {}).get(old_name_clean, new_name) != new_name:
                    print(f"Warning: {old_name_clean} for organism type {organism_type} is associated with multiple new names: {ribo_name_organism_conversions.get(organism_type, {}).get(old_name_clean)} and {new_name}. Keeping the last one.")
                ribo_name_organism_conversions.setdefault(organism_type, {})[old_name_clean] = new_name
# ensure preferred conversions are included (and take precedence over non-preferred ones in the regex)
ribo_name_conversions.update(ribo_name_conversions_preferred)

ribo_subunit_regex = r'[abmu][SL]\d+m{0,1}\w{0,1}|M{0,1}RP[SLP]\d+' # SC naming scheme, plus scheme from https://bangroup.ethz.ch/research/nomenclature-of-ribosomal-proteins.html
# update regex to also include each entry in set(ribo_name_conversions.keys()) | set(ribo_name_conversions.values()) as separate words
conversion_names_regex = r'|'.join(sorted([re.escape(name) for name in set(ribo_name_conversions.keys()) | set(ribo_name_conversions.values()) if not re.match(ribo_subunit_regex, name)], key=lambda x: len(x), reverse=True)) # sort by length in reverse order to avoid partial matches (e.g., RPS6A matching RPS6)
ribo_subunit_regex = r'\b(' + conversion_names_regex + '|' + ribo_subunit_regex + r')\b' if conversion_names_regex != '' else ribo_subunit_regex
for i in fluxes:
    r,vtype,val_raw = i.split('\t')
    val = float(val_raw) / nscale2
    fluxes_new.append('\t'.join([r, vtype, str(val)]))
    rxns_used.append("'" + r + "'")
    tag,rxn_base_id,rxn_dir,enz_id = extract_details_from_rxnid(r)
    # if rxn (as determined by characters after last "-") isn't SPONT or UNKNOWN
    if enz_id not in spont_rxn_suffixes and tag == 'RXN':
        for enzload in rxn_to_enzloads.get(r, []):
            enz = enzload.split('-', 1)[1].rsplit('-', 1)[1]
            enz_rxn_dict.setdefault(enz, []).append(enzload)
            enzyme_rxns.append("'" + enzload + "'")
    elif tag in ['PROSYN','PTM']:
        pro = r.split('-',1)[1]
        # find df_prot matches by id, entry, or uniprot
        df_prot_match = df_prot[(df_prot['id'] == pro) | (df_prot['gene_src'] == pro) | (df_prot['Entry'] == pro) | (df_prot['uniprot'] == pro)]
        if len(df_prot_match) == 0:
            print('No match found in df_prot for', pro)
        if tag in ['PROSYN','PTM']:
            if tag == 'PROSYN':
                waste_rxn = r.replace('PROSYN-','PROWASTE-',1)
                # calculate amino acid composition, for use in the GSM model (if different/needed)
                # find the amino acid sequence of the protein
                aa_seq = df_prot.loc[df_prot_match.index[0], 'sequence']
                for aa in aa_seq:
                    # check if it's a letter
                    if aa.isalpha():
                        # add val to aa_comp, since for every mol of protein made, 1 mol of each aa in its sequence is made
                        if aa in aa_comp.keys():
                            # aa_comp[aa]['total flux'] += val_raw
                            aa_comp[aa] += float(val_raw)
                        else:
                            # aa_comp[aa] = {'total flux':val_raw}
                            aa_comp[aa] = float(val_raw)
                # enter data into kribo calculations
                # do this separately for each ribosome type and establish a min and max range for any proteins that can be translated by multiple ribosome types, to be used in kribo calculation
                if len(pro_trans_dict[pro]) > 1:
                    for ribosome in pro_trans_dict[pro]:
                        V_ribo_K_ribo[ribosome]['max'] += mu * protein_length[pro] * float(val)
                else:
                    for value in ['min','max']:
                        V_ribo_K_ribo[pro_trans_dict[pro][0]][value] += mu * protein_length[pro] * float(val)
            elif tag == 'PTM':
                # check if its val equals that of the flux_dict entry for its associated PTMWASTE reaction, to see if it's not used by any enzymes
                waste_rxn = r.replace('PTM-','PTMWASTE-',1)
            if val != 0:
                protein_id = met_ID.from_name(r.replace('PTM-','PTMPRO-')).base_id
                if 'Gene Ontology IDs' in df_prot.columns and pd.notna(df_prot.loc[df_prot_match.index[0], 'Gene Ontology IDs']):
                    go_terms = set(df_prot.loc[df_prot_match.index[0], 'Gene Ontology IDs'].strip().split('; '))
                    keywords = set(df_prot.loc[df_prot_match.index[0], 'Keywords'].strip().split(';')) if 'Keywords' in df_prot.columns and pd.notna(df_prot.loc[df_prot_match.index[0], 'Keywords']) else set()
                    subcellular_locations = set([loc for loc in re.sub(r'; .*?(\.|$)','. ',re.sub(r'\s+', ' ', re.sub(r' \{.*?\}|\.$', '', df_prot.loc[df_prot_match.index[0], 'Subcellular location [CC]'].replace('SUBCELLULAR LOCATION: ','')))).strip().split('. ')]) if 'Subcellular location [CC]' in df_prot.columns and pd.notna(df_prot.loc[df_prot_match.index[0], 'Subcellular location [CC]']) else set()
                    gene_names = set(df_prot.loc[df_prot_match.index[0], 'Gene Names'].strip().split(' ')) if 'Gene Names' in df_prot.columns and pd.notna(df_prot.loc[df_prot_match.index[0], 'Gene Names']) else set()
                    potential_ribo_subunit_names = set()
                    subunit_groups_initial = set()
                    for col in ['Protein families', 'Protein names', 'Gene Names']:
                        colval = df_prot.get(col, {}).get(df_prot_match.index[0], '')
                        colval = colval if isinstance(colval, str) else ''
                        subunit_ids = set(re.findall(r'\b(\d+S)\b[^\t]*(?i:ribo)', colval))
                        if len(subunit_ids) > 1:
                            print(f"Warning: multiple potential ribosome subunit IDs found in {col} for {protein_id}: {subunit_ids}. This may be due to multiple subunits being listed or due to parsing errors. Please check the formatting of the {col} column for this protein and consider manually assigning ribosome subunit membership if needed.")
                        if 'ribo' not in colval.lower() and col != 'Gene Names':
                            continue
                        for match in re.findall(ribo_subunit_regex, colval):
                            # convert name using subunit_ids if possible (e.g., if subunit_id is 40S or 60S, use ribo_name_organism_conversions entries except for bacteria; if it's 37S or 54S, use yeast ribosome-specific entries; it it's 28S, use mammal entries; otherwise, leave it alone
                            if subunit_ids:
                                converted_names = set()
                                subunit_id = list(subunit_ids)[0]
                                if subunit_id in ['40S', '60S']:
                                    if RBA_defaults.organism_type.lower() not in organism_types.difference(['bacteria', 'archaea']):
                                        for organism in ribo_name_organism_conversions.keys():
                                            if organism in ['bacteria', 'archaea']:
                                                continue
                                            if match in ribo_name_organism_conversions.get(organism, {}).keys():
                                                converted_names.add(ribo_name_organism_conversions.get(organism, {}).get(match, match))
                                                break
                                    else: # assume if it's eukaryotic, that the entries are formatted for the organism of interest (important since)
                                        # print(f"Subunit ID {subunit_id} found for {match} in {col} for {protein_id}, but organism type is {RBA_defaults.organism_type} and ribo_name_organism_conversions entries for this subunit ID are not formatted for this organism type. Adding unconverted name as potential ribosome subunit name, but please check the formatting of the {col} column for this protein and consider manually assigning ribosome subunit membership if needed.")
                                        converted_names.add(ribo_name_organism_conversions.get(RBA_defaults.organism_type.lower(), {}).get(match, match))
                                elif subunit_id in ['37S', '54S']:
                                    # print(f"Subunit ID {subunit_id} found for {match} in {col} for {protein_id}, but since this subunit ID is specific to yeast ribosomes, converting name using yeast-specific entries in ribo_name_organism_conversions. Please check the formatting of the {col} column for this protein and consider manually assigning ribosome subunit membership if needed.")
                                    converted_names.add(ribo_name_organism_conversions.get('yeast', {}).get(match, match))
                                elif subunit_id == '28S':
                                    # print(f"Subunit ID {subunit_id} found for {match} in {col} for {protein_id}, but since this subunit ID is specific to mammalian ribosomes, converting name using mammal-specific entries in ribo_name_organism_conversions. Please check the formatting of the {col} column for this protein and consider manually assigning ribosome subunit membership if needed.")
                                    converted_names.add(ribo_name_organism_conversions.get('mammal', {}).get(match, match))
                                    # print(f"Converted name for {match} in {col} for {protein_id} based on subunit ID {subunit_id}: {converted_names}. Please check the formatting of the {col} column for this protein and consider manually assigning ribosome subunit membership if needed.")
                                else:
                                    # print(f"Subunit ID {subunit_id} found for {match} in {col} for {protein_id}, but this subunit ID is not recognized. Adding unconverted name as potential ribosome subunit name, but please check the formatting of the {col} column for this protein and consider manually assigning ribosome subunit membership if needed.")
                                    converted_names.add(match)

                                if len(converted_names) > 1:
                                    print(f"Warning: multiple converted names found for {match} with subunit ID {subunit_id} for {protein_id}: {converted_names}. This may be due to multiple matches in the ribo_name_organism_conversions dictionary. Please check the entries in ribo_name_organism_conversions for this subunit ID and consider manually assigning ribosome subunit membership if needed.")
                                potential_ribo_subunit_names.update(converted_names)
                            else:
                                potential_ribo_subunit_names.add(match)
                            
                    if 'GO:0003735' in go_terms or 'Ribosomal protein' in keywords or len(potential_ribo_subunit_names) > 0:
                        for ribosome in V_ribo_K_ribo.keys():
                            ribosome_confirmed = False
                            subunit_groups = set()  # Initialize as a list
                            subunit_groups.add('SSU') if 'GO:0015935' in go_terms else None
                            subunit_groups.add('LSU') if 'GO:0015934' in go_terms else None
                            if ribosome in RBA_defaults.mito_compartments:
                                ssu_terms = set(['GO:0005763'])
                                lsu_terms = set(['GO:0005762'])
                                # check subcellular location for mitochondrion
                                if go_terms.intersection(['GO:0005761']) or 'Mitochondrion' in subcellular_locations:
                                    ribosome_confirmed = True
                            elif ribosome == RBA_defaults.cytosol:
                                ssu_terms = set(['GO:0022627', 'GO:0030873', 'GO:0030499', 'GO:0009283', 'GO:0005843'])
                                lsu_terms = set(['GO:0022625', 'GO:0030498', 'GO:0030872', 'GO:0009282', 'GO:0005842'])
                                if go_terms.intersection(['GO:0022626']) or 'Cytoplasm' in subcellular_locations or RBA_defaults.cytosol in df_prot.loc[df_prot_match.index[0], 'subloc_assigned'].split(','):
                                    ribosome_confirmed = True
                            else:
                                ssu_terms = set(['GO:0000314'])
                                lsu_terms = set(['GO:0000315'])
                            if go_terms.intersection(lsu_terms):
                                subunit_groups.add('LSU')
                                ribosome_confirmed = True
                            elif go_terms.intersection(ssu_terms):
                                subunit_groups.add('SSU')
                                ribosome_confirmed = True
                            if ribosome_confirmed:
                                subunit_groups = list(subunit_groups) if subunit_groups else ['unknown']
                                for subunit_group in subunit_groups:
                                    ribo_subunit_candidates[ribosome][subunit_group][protein_id] = {'flux': val, 'GO_terms': sorted(list(go_terms))}
                                    if len(potential_ribo_subunit_names) > 0:
                                        potential_ribo_subunit_names_converted = set()
                                        for name in potential_ribo_subunit_names:
                                            # name = interpret_subunit_name(name, return_format='string, paralog-agnostic') # convert to a standard format that ignores paralog information, since it's often not clear from the proteomics data which paralog is which and they may be interchangeable
                                            if name in ribo_name_conversions.keys():
                                                potential_ribo_subunit_names_converted.add(ribo_name_conversions[name])
                                            elif name in ribo_name_conversions.values():
                                                potential_ribo_subunit_names_converted.add(name)
                                        potential_ribo_subunit_names_converted = sorted(list(potential_ribo_subunit_names_converted))
                                        ribo_subunit_candidates[ribosome][subunit_group][protein_id]['potential_ribo_subunit_names'] = potential_ribo_subunit_names_converted
                                        for subunit_name in potential_ribo_subunit_names_converted:
                                            subunit_name_converted = interpret_subunit_name(subunit_name, return_format='string, paralog-agnostic')
                                            # subunit_dict = interpret_subunit_name(subunit_name)
                                            # filter out subunit names that don't match the assigned subunit groups, since there's less evidence for it
                                            # if subunit_group != 'unknown' and subunit_dict['subunit_group'] != subunit_group[0]:
                                            #     continue
                                            ribo_subunit_candidates[ribosome]['composition_dict'].setdefault(subunit_name_converted, []).append(protein_id)

                                    if flux_dict.get(waste_rxn) and flux_dict.get(waste_rxn, 0) != val:
                                        ribo_subunit_candidates[ribosome][subunit_group][protein_id]['disclaimers'] = ["Protein used for other purposes, making it hard to determine how much is allocated to ribosome (if any)."]
            # Optional per-protein amino acid composition calculations
            # json.dump(aa_comp, open(os.path.join(path_out, 'enz_alloc_aa_comp_'+enz+'.json'), 'w'), indent=4)

for ribo in ribo_subunit_candidates.keys():
    min_flux_needed = V_ribo_K_ribo[ribo]['min'] / RBA_defaults.k_ribo_max
    for subunit_group in ribo_subunit_candidates[ribo].keys():
        if subunit_group in ['composition_dict']:
            continue
        # if flux * RBA_defaults.k_ribo_max < V_ribo_K_ribo for this ribosome type, list such proteins as potential paralogs since there can't be 1 present per ribosome
        for pro, details in ribo_subunit_candidates[ribo][subunit_group].items():
            ribo_subunit_candidates[ribo][subunit_group][pro]['max copies per ribosome'] = f"{(details['flux']/min_flux_needed):.2f}"
            if details['flux'] < min_flux_needed:
                ribo_subunit_candidates[ribo][subunit_group][pro].setdefault('disclaimers', []).append(f"Likely interchangeable with other proteins in the same subunit group, since no more than {(details['flux']/min_flux_needed):.2f} molecules per ribosome can be provided.")
    # in the composition dict, convert all values to tuples; calculate total flux for each subunit based on the sum of fluxes of its candidate proteins
    for subunit, proteins in ribo_subunit_candidates[ribo]['composition_dict'].copy().items():
        if type(subunit) != tuple:
            new_key = tuple(dict.fromkeys(proteins))
            expected_copies_per_ribosome = 2 if subunit in ['P1','P2'] else 1 # ribosomal stalk protein composition based on https://doi.org/10.1093/nar/gkt636 and https://doi.org/10.1093/nar/gks036
            ribo_subunit_candidates[ribo]['composition_dict'][new_key] = {"expected_copies_per_ribosome": expected_copies_per_ribosome, 
                                                                          "total_flux": sum(ribo_subunit_candidates[ribo][subunit_group].get(pro, {}).get('flux', 0) for pro in proteins for subunit_group in ['SSU','LSU','unknown'])}
            # check if total flux surpasses min_flux_needed for this ribosome type, and if not, remove it
            if ribo_subunit_candidates[ribo]['composition_dict'][new_key]['total_flux'] < min_flux_needed * expected_copies_per_ribosome:
                ribo_subunit_candidates[ribo]['composition_dict'].pop(new_key)
            ribo_subunit_candidates[ribo]['composition_dict'].pop(subunit)

# for each ribosome, find the subunit in composition_dict with the lowest total flux, and calculate kribo based on that flux and the assumption that all of it is allocated to the ribosome (i.e. it's not a paralog with other proteins in the same subunit group that could be interchangeable)
for ribo in ribo_subunit_candidates.keys():
    min_flux = float('inf')
    for subunit, details in ribo_subunit_candidates[ribo]['composition_dict'].items():
        if details['total_flux'] / details['expected_copies_per_ribosome'] < min_flux:
            min_flux = details['total_flux'] / details['expected_copies_per_ribosome']
    # update composition_dict values by dividing by min_flux, to get relative abundances of each subunit compared to the one with lowest flux (round down to nearest integer since you can't have a fraction of a protein subunit in a ribosome)
    for subunit, details in ribo_subunit_candidates[ribo]['composition_dict'].items():
        ribo_subunit_candidates[ribo]['composition_dict'][subunit] = {'max copies per ribosome': int(details['total_flux'] / min_flux), 'expected copies per ribosome': details['expected_copies_per_ribosome']}
    kribo = V_ribo_K_ribo[ribo]['max'] / min_flux if min_flux > 0 else 0
    ribo_subunit_candidates[ribo]['kribo'] = float(kribo)

# divide aa_comp values by biomass flux to get new GSM coefficients, in case new ones are needed
for aa in aa_comp.keys():
    aa_comp[aa] = aa_comp[aa]/nscale2
    aa_biom_coeff = aa_comp[aa]/(mu)
    # find matching tRNA name in aa_map, and create draft of biomass rxn from it (in case it differs from GSM)
    biomass_from_RBA[aa_map[aa_map['aa_abbv'] == aa.upper()]['tRNA_in'].to_list()[0].split('MET-')[1]] = -aa_biom_coeff
    biomass_from_RBA[aa_map[aa_map['aa_abbv'] == aa.upper()]['tRNA_out'].to_list()[0].split('MET-')[1]] = aa_biom_coeff
    
# calculate mol % of each aa in aa_comp
# for aa in aa_comp.keys():
#     aa_comp[aa]['mol %'] = aa_comp[aa]['total flux']/sum([aa_comp[aa]['total flux'] for aa in aa_comp.keys()])
# store aa composition, biomass details in a file
json.dump(biomass_from_RBA, open(os.path.join(path_out, 'enz_alloc_biomass_from_RBA_partial.json'), 'w'), indent=-1, sort_keys=True)
json.dump(aa_comp, open(os.path.join(path_out, 'enz_alloc_total_aa_fluxes_for_prot_unscaled.json'), 'w'), indent=-1, sort_keys=True)
# find lowest flux out of all vtrans values in df_data and all fluxes in fluxes_new, as well as all predicted ENZLOAD fluxes based on EnzLoadConstraint values
all_nonzero_fluxes = df_data[df_data['vtrans (mmol/gDW/h)'] > 0]['vtrans (mmol/gDW/h)'].to_list()
# all_nonzero_fluxes = df_data[df_data['vtrans (mmol/gDW/h)'] > 0]['vtrans (mmol/gDW/h)'].to_list() + [float(i.split('\t')[2]) for i in fluxes_new if float(i.split('\t')[2]) > 0] \
# + [float(i.split('\t')[2])/len(enz_rxn_dict[i.split('\t')[2]]) for i in fluxes_new if str(i.split('\t')[0]).startswith('ENZSYN-')]
enzloads_predicted = dict()
for i in fluxes_new:
    r,vtype,val = i.split('\t')
    val = float(val) * nscale2
    prefix,enz_draft = r.split('-')[0],r.split('-')[-1]
    if val > 0:
        all_nonzero_fluxes.append(val)
        if enz_draft not in spont_rxn_suffixes:
            if prefix == 'RXN':
                for enzload in rxn_to_enzloads.get(r, []):
                    enzloads_predicted[enzload] = 0
            elif prefix == 'ENZSYN':
                # find 'enz' value in sij_dict for this rxn, and use it to predict the flux for its associated ENZLOAD reaction, which is a proxy for enzyme abundance; add predicted ENZLOAD flux to list of all_nonzero_fluxes to determine min and max flux values for solver
                enzs = [enz for enz in sij_dict[r].keys() if enz.startswith('ENZ-')]
                if len(enzs) > 1:
                    print('Multiple enzymes found in sij_dict for', r, ':', enzs)
                enz_full_id = enzs[0]
                enz = met_ID.from_name(enz_full_id).base_id
                if enz in enz_rxn_dict.keys():
                    v_enzload_pred = val/len(enz_rxn_dict[enz])
                    # add predicted ENZLOAD flux, in case the solver treated it as 0
                    if v_enzload_pred > 0:
                        all_nonzero_fluxes.append(v_enzload_pred)
                    else:
                        print('Predicted ENZLOAD flux for each rxn using', enz, 'is 0')
                else:
                    rxns_using_enz = [r for r in sij_dict.keys() if r.startswith('ENZLOAD-') and enz_full_id in sij_dict[r].keys()]
                    if len(rxns_using_enz) > 0:
                        print(f"Enzyme '{enz}' used in the following rxns but not in metabolic rxns: {rxns_using_enz}")
                    else:
                        print(f"Enzyme '{enz}' not used in any ENZLOAD rxns")
            elif prefix == 'ENZLOAD':
                enzloads_predicted[r] = val
# check if any ENZLOAD fluxes were predicted to be 0 despite their associated rxn fluxes having nonzero values
for enzload,v in enzloads_predicted.items():
    if v == 0:
        print('Predicted ENZLOAD flux for', enzload, 'is 0 but its associated rxn has nonzero flux')
# all_nonzero_fluxes = list(df_data['vtrans (mmol/gDW/h)']) + [float(i.split('\t')[2]) for i in fluxes_new] \
    # + [float(i.split('\t')[2])/len(enz_rxn_dict[i.split('\t')[2]]) for i in fluxes_new if str(i.split('\t')[0]).startswith('ENZSYN-')]
v_min = min(all_nonzero_fluxes)
v_max = max(all_nonzero_fluxes)
print('v_min:', v_min, 'v_max:', v_max)
print('Make sure both values fit within the solver tolerance range')

fname = os.path.join(path_out, 'enz_alloc.flux.txt')
with open(fname, 'w') as f:
    f.write('\n'.join(fluxes_new))

# output composition dict for each ribosome to YAML file for use in kribo calculation
ribo_composition_file = 'ribosome_composition.yaml'
ribo_composition_path = os.path.join(path_out, ribo_composition_file)
with open(ribo_composition_path, 'w') as f:
    yaml.dump({ribo:{'k_ribo': ribo_subunit_candidates[ribo]['kribo'], 'subunits': {' OR '.join(proteins): details["expected copies per ribosome"] for proteins, details in ribo_subunit_candidates[ribo]['composition_dict'].items()}} for ribo in ribo_subunit_candidates.keys()}, f, default_flow_style=False)
if not any(os.path.isfile(v) and os.path.exists(f"{model_root_path}build_model/{v}") for k, v in getattr(RBA_defaults, 'ribo_dict', {}).items()):
    if not os.path.isfile(f"{model_root_path}build_model/input/{ribo_composition_file}"):
        os.system(f"cp {ribo_composition_path} {model_root_path}build_model/input/{ribo_composition_file}")
for ribo in ribo_subunit_candidates.keys():
    # convert composition dict keys from tuples to strings for better readability in JSON output
    ribo_subunit_candidates[ribo]['composition_dict'] = {str(proteins): details for proteins, details in ribo_subunit_candidates[ribo]['composition_dict'].items()}
# output ribo_subunit_candidates to file for use in kribo calculation
json.dump(ribo_subunit_candidates, open(os.path.join(path_out, 'enz_alloc_ribo_subunit_candidates.json'), 'w'), indent=2, sort_keys=True)

used_subunits = set([p for ribo in ribo_subunit_candidates.keys() for proteins in ribo_subunit_candidates[ribo]['composition_dict'].keys() for p in eval(proteins)])
# output candidates that weren't used in the final ribosome composition dict to a separate file for reference
unused_candidates_per_subunit = {ribo:
                     {subunit_group:
                      {pro: details for pro, details in ribo_subunit_candidates[ribo][subunit_group].items() 
                       if pro not in [p for proteins in ribo_subunit_candidates[ribo]['composition_dict'].keys() for p in eval(proteins)]} 
                      for subunit_group in ribo_subunit_candidates[ribo].keys() if subunit_group not in {'composition_dict', 'kribo'}} 
                    for ribo in ribo_subunit_candidates.keys()}
candidates_never_used = {ribo:
                            {subunit_group:
                            {pro: details for pro, details in ribo_subunit_candidates[ribo][subunit_group].items() 
                            if pro not in used_subunits}
                            for subunit_group in ribo_subunit_candidates[ribo].keys() if subunit_group not in {'composition_dict', 'kribo'}} for ribo in ribo_subunit_candidates.keys()}
json.dump(unused_candidates_per_subunit, open(os.path.join(path_out, 'enz_alloc_ribo_subunit_candidates_unused_per_subunit.json'), 'w'), indent=2, sort_keys=True)
json.dump(candidates_never_used, open(os.path.join(path_out, 'enz_alloc_ribo_subunit_candidates_never_used.json'), 'w'), indent=2, sort_keys=True)

for ribo in ribo_subunit_candidates.keys():
    print(ribo, [p for proteins in ribo_subunit_candidates[ribo]['composition_dict'].keys() for p in eval(proteins)])

print(f'V_ribo * K_ribo for each ribosome type: {V_ribo_K_ribo}')
# run kribo_bisection.gms
cmds = ['cd ' + path_out,
        'module load gams',
        f"gams kribo_bisection.gms --nscale={str(nscale)} --kapp_max={kapp_max} {output_redirect_str}"]
#os.system('\n'.join(cmds))
#stop_if_run_failed('./min_flux_violation/kribo_bisection.modelStat.txt')

with open('./min_flux_violation/enz_alloc.flux_gamsscaled.txt') as f:
    rxns_used = f.read().split('\n')
# find flux of total protein going into biomass: BIOSYN-PROTTOBIO
# make dict of all PROWASTE rxns (except PROWASTE-TOTALPROTEIN) with name as key and flux as value
prowaste_used = dict()
prosyn_used = dict()
rxn_ids_used = set([i.split('\t')[0].replace("'",'') for i in rxns_used])
ribo_exists = any(i.startswith('RIBOSUBSYN-') for i in rxn_ids_used)
for i in rxns_used:
    i_twocol = i.replace('	v	 ','\t').replace('\tv\t','\t')
    if len(i_twocol.split('\t')) < 2:
        print(i_twocol)
        continue
    rxn,v = i_twocol.split('\t')[:2]
    v = float(v)/nscale2
    tag,rxn_base_id,direc,enz = extract_details_from_rxnid(rxn)
    i_twocol = rxn+'\t'+str(v)
    waste_syn = {'PROWASTE': 'PROSYN', 'PTMWASTE': 'PTM'}
    if tag in waste_syn.keys() and rxn != 'PROWASTE-TOTALPROTEIN':
        syn_id = rxn.replace(tag, waste_syn.get(tag, ''), 1)
        prowaste_used[syn_id] = prowaste_used.get(syn_id, 0) + v
    elif tag == 'PROSYN':
        prosyn_used[rxn] = v 
    elif rxn == 'BIOSYN-PROTTOBIO':
        v_prot_to_bio = float(i_twocol.split('\t')[1]) * nscale2 # nscale2 added back for consistency with mass fractions

# Add new set of constraints to reflect unused but observed protein production, similar to dummy protein (for every prot, v('PROSYN-prot') >= prowaste_used('PROWASTE-prot') / v_prot_to_bio * %nscale%)

## UNUSED: make separate file for defaults, so not every options file needs to be updated when the model changes
## pro_measured_unused_default_path = '../../../../prosyn_constraints_measured_but_unused_default.txt'
## if os.path.exists(pro_measured_unused_default_path):

pro_measured_unused_path = 'prowaste_from_pro_with_data.txt'
pro_predicted_unused_path = 'prowaste_from_pro_without_data.txt'
pro_measured_unused_path_without_ribo = 'prowaste_from_pro_with_data_without_ribo.txt'
pro_predicted_unused_path_without_ribo = 'prowaste_from_pro_without_data_without_ribo.txt'

# predicted flux allocated to ribosome subunits for each protein
ribo_flux_by_protein = dict()
for ribo_id, ribo_info in ribo_subunit_candidates.items():
    for subunit_group, details in ribo_info['composition_dict'].items():
        for protein in eval(subunit_group):
            for subunit in ['SSU', 'LSU', 'unknown']:
                if protein in ribo_info[subunit]:
                    ribo_flux_by_protein[protein] = min(float(ribo_info[subunit][protein].get('flux', 0)), (V_ribo_K_ribo[ribo_id]['min'] / ribo_info['kribo']) * details.get('expected copies per ribosome', 1)) # cap the flux allocated to the ribosome for this protein at the minimum flux needed for the ribosome divided by kribo, multiplied by the expected number of copies of this protein per ribosome, since it's unlikely that more than that would be allocated to the ribosome; this is a conservative approach that may underestimate the flux allocated to the ribosome for some proteins but ensures we don't overestimate it for any proteins

with open(pro_measured_unused_path,'w') as f, open(pro_predicted_unused_path,'w') as f2, open(pro_measured_unused_path_without_ribo,'w') as f4, open(pro_predicted_unused_path_without_ribo,'w') as f5, open('prot_mass_flux.txt','w') as f3:
    f3.write('/'+str(v_prot_to_bio)+'/')
    # for each protein in PROWASTE-used, allow flux through PROWASTE-PROT to be as high as infinity (to avoid scaling issues)
    # Also, make the proteome mass fraction comprised of each protein at least equal to the mass fraction comprised of unused copies of it 
    ## (i.e., PROSYN * MW >= ptot * PROWASTE_0 * MW / ptot_0)
    ### divide both sides by MW: PROSYN >= ptot * PROWASTE_0 / ptot_0)
    for prot in prowaste_used.keys():
        prot_info = met_ID.from_name(prot)
        prot_id = prot_info.base_id + '_' + prot_info.compartment if prot_info.compartment else prot_info.base_id
        if prot not in rxn_ids_used:
            print(f"Warning: rxn ID '{prot}' not found in model reactions. Skipping constraint for this protein. Please check the formatting of the PROWASTE reaction IDs and ensure they match the corresponding PROSYN reaction IDs for accurate constraint generation.")
            continue
        constraint = "'"+prot+"' "+str(prowaste_used[prot])+"\n"
        # Only subtract predicted ribosomal contributions if the ribosomes aren't already defined, to avoid subtracting them twice
        without_ribo_val = max(prowaste_used[prot] - ribo_flux_by_protein.get(prot_info.base_id, 0), 0) if not ribo_exists else prowaste_used[prot]
        without_ribo_val = str(without_ribo_val) if without_ribo_val == prowaste_used[prot] else f"{without_ribo_val} /* subunit */"
        without_ribo_constraint = "'"+prot+"' "+str(without_ribo_val)+"\n"
        # if protein is in the proteome data (path_data), add to f2
        if prot_id in df_data_full['uniprot'].values.tolist() + df_data_full.index.tolist():
            f.write(constraint)
            f4.write(without_ribo_constraint)
        else:
            f2.write(constraint)
            f5.write(without_ribo_constraint)
