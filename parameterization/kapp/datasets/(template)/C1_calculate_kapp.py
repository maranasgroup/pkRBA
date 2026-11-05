# update model-specific settings in kapp_options.py
load_model = False
from kapp_options import *
min_flux_cutoff = 0
# Change to script directory to use relative paths
os.chdir(os.path.dirname(os.path.abspath(__file__))) 
# tolerance value added to assist with rounding kapps; set to 0 if not needed
tol = 0 
round_up_if_infeas = False
round_up_default = False
if len(sys.argv) > 1:
    round_up_default = bool(int(sys.argv[1]))

def kapp_conversion(value, round_up=round_up_default):
    """kapp processing,optionally rounding its final digit up."""
    text = str(value)
    # round up last digit to prevent rounding issues with very low kapps; e.g., if kapp is 1e-5, round to 1.1e-5
    if round_up and '.' in text:
        value = float(text[:-1] + str(int(text[-1]) + 1))
    return value

def read_lines(path, replace_char="'", skip='/', strip=True):
    with open(path) as stream:
        lines = [line.replace(replace_char, '') for line in stream.read().splitlines()]
    if strip: lines = [line.strip() for line in lines]
    lines = [line for line in lines if line not in (skip, '')]
    return lines

# remove location to see if enzyme is identical to another
def write_lines(path, lines):
    with open(path, 'w') as stream:
        stream.write('\n'.join(lines))

def direction(flux, label):
    if flux == 0:
        print('rval == 0, check enzyme ' + label)
    return 'FWD' if flux > 0 else 'REV'

def enz_load(enz, raw_flux, sij):
    return sum(value for rxn, value in raw_flux.items()
               if 'ENZ-' + enz in sij.get(rxn, {}))

def rid(rxn, flux, enz):
    return f"ENZLOAD-{rxn}_{direction(flux, enz + ' and ' + rxn)}-{enz}"

def enzload_id(rxn):
    return rxn.replace('RXN-', 'ENZLOAD-', 1)

def save_csv(path, values):
    with open(path, 'w') as stream:
        stream.write('rxnid,median\n')
        for key, vals in sorted(values.items()):
            stream.write(','.join(map(str, [key, min(vals), max(vals),
                                            np.median(vals), vals])) + '\n')

# Read model and flux inputs.
sij = {}
for line in read_lines(sij_path, replace_char=""):
    met, rxn, value = line.replace("'.'", "' '").split()
    sij.setdefault(rxn[1:-1], {})[met[1:-1]] = float(value)
df_enz = read_spreadsheet('../../../../build_model/model/ENZYME_stoich_curation.tsv')
# set4 handles cases of multiple enzymes corresponding to multiple rxns (many-to-many)
# columns: 
#   enzyme (for all enz that correspond to a rxn)
#   rxn: 
#       disableManually for any I don't want to count
#       noRxnFlux_resolveByManualCheck: behaves identically to disableManually, but perhaps serves as a warning
#   rule: 
#       divmaxenz for cases where I want to treat the [enzyme] as equaling its maximum value
set4_path = './kapp_ambiguousLoad_case_resolve_common.txt'
res_metab = RBA_result(biom_id=biom_id)
flux_path = './min_flux_violation/enz_alloc.flux_gamsscaled.txt'
# check if any enzymatic rxns with no proteomics data were used
rxns_essential_NP = [line.split('\t')[0] for line in
                     open('min_flux_violation/enz_alloc.flux_essential_with_no_prodata_unscaled.txt')
                     .read().splitlines() if line]

prowaste_used, prosyn_used, enzsyn_used, enz_rxns_without_enzload = {}, {}, [], {}
v_prot_to_bio = 0
for line in open(flux_path).read().splitlines():
    parts = line.replace('\tv\t', '\t').replace('\tv\t', '\t').split('\t')
    if len(parts) < 2:
        print(line); continue
    rxn, value = parts[:2]; value = float(value) / nscale2
    tag, _, _, enz = extract_details_from_rxnid(rxn)
    if tag == 'PROWASTE' and rxn != 'PROWASTE-TOTALPROTEIN': prowaste_used[rxn] = value
    elif tag == 'PROSYN': prosyn_used[rxn] = value
    elif tag == 'ENZSYN': enzsyn_used.append(f'{rxn}\t{value}')
    elif rxn == 'BIOSYN-PROTTOBIO': v_prot_to_bio = value
    elif tag == 'RXN' and enz not in spont_rxn_suffixes: enz_rxns_without_enzload[rxn] = value

write_lines('./kapp_input_enzsyn.txt', enzsyn_used)
res_metab.load_raw_flux(flux_path, nscale=nscale2)
res_metab.calculate_metabolic_flux()
res_metab.make_escher_csv('./flux_metab.escher.csv')
res_esyn = RBA_result(biom_id=biom_id, twocol_format=True)
res_esyn.load_raw_flux('./kapp_input_enzsyn.txt')
all_metab_rxns, all_GSM_rxns = read_lines(metab_rxns_path), read_lines(path_model + 'SM_rxn_ids.txt')
all_rxn_enzload_pairs = dict()
for r in sij:
    if not r.startswith('ENZLOAD-'): continue
    tag, base, dir, enz = extract_details_from_rxnid(r)
    all_rxn_enzload_pairs.setdefault(f"RXN-{base}_{dir}", []).append(r)
mu = res_metab.growth_rate
print('Growth rate:', mu)

# Build enzyme/reaction relationships and classify their load-sharing patterns.
rxndict = {rxn: [] for rxn, flux in res_metab.metabolic_flux.items() if abs(flux) > 0}
enzdict = {}
for rxn, value in res_esyn.raw_flux.items():
    if value > 0:
        enz = met_ID.from_name(next(e for e in sij[rxn] if e.startswith('ENZ-')))
        key = str(met_ID(base_id=enz.base_id, compartment=enz.compartment)).removeprefix('MET-')
        enzdict[key] = []
for row in df_enz.itertuples():
    if row.rxn_src not in rxndict: continue
    if row.enz in spont_rxn_suffixes: rxndict[row.rxn_src].append('zeroCost')
    elif row.enz in enzdict:
        rxndict[row.rxn_src].append(row.enz); enzdict[row.enz].append(row.rxn_src)
rxndict = {k: set(v) for k, v in rxndict.items()}
zero_cost = {k for k, v in rxndict.items() if v == {'zeroCost'}}
rxndict = {k: v for k, v in rxndict.items() if k not in zero_cost and 'zeroCost' not in v}
enzdict = {k: set(v) for k, v in enzdict.items() if v and k not in zero_cost}
shared = set().union(*(v for v in rxndict.values() if len(v) > 1))
multi = {k for k, v in enzdict.items() if len(v) > 1}
sets = [set(enzdict) - shared - multi, multi - shared, shared - multi, shared & multi]
for number, values in enumerate(sets, 1):
    write_lines(f'./kapp_set{number}_' +
                ['one_to_one.txt', 'one_to_many.txt', 'many_to_one.txt', 'many_to_many.txt'][number-1],
                sorted(values))
set1, set2, set3, set4 = sets
set4 = set4 - set(rxns_to_ignore_for_kapps)
with open(set4_path, 'w') as stream:
    stream.write('Enzyme\tReaction\tRule\n')
    for enz in sorted(set4): stream.write(f'{enz}\t{",".join(enzdict[enz])}\t\n')

# Calculate the one-to-one, one-to-many, and many-to-one cases.
kapp, kapps_rba = {}, {}
for enz in set1 | set2:
    rxns = list(enzdict[enz]); total = sum(abs(res_metab.metabolic_flux[r]) for r in rxns)
    load = enz_load(enz, res_esyn.raw_flux, sij)
    for rxn in rxns:
        key = rid(rxn, res_metab.metabolic_flux[rxn], enz)
        value = mu * (abs(res_metab.metabolic_flux[rxn]) if enz in set1 else total) / load
        kapp[key], kapps_rba[key] = value / 3600, kapp_conversion(value)
for rxn in set().union(*(enzdict[e] for e in set3)):
    flux = res_metab.metabolic_flux[rxn]
    load = sum(enz_load(e, res_esyn.raw_flux, sij) for e in rxndict[rxn] if e != 'zeroCost')
    for enz in rxndict[rxn] - {'zeroCost'}:
        key, value = rid(rxn, flux, enz), mu * abs(flux) / load
        kapp[key], kapps_rba[key] = value / 3600, kapp_conversion(value)

# Manual many-to-many resolutions.
mapper = {}
for line in read_lines(set4_path, strip=False)[1:]:
    enz, rxns, rule = line.split('\t')
    if rxns not in ('noRxnFlux_resolveByManualCheck', 'disableManually'):
        mapper[enz] = (rxns, rule)
for enztext, (rxntext, rule) in mapper.items():
    enzs = enztext.split(',')
    measured = [e for e in enzs if enz_load(e, res_esyn.raw_flux, sij)]
    if rule == 'divmaxenz':
        load = max(enz_load(e, res_esyn.raw_flux, sij) for e in measured)
    else: 
        load = sum(enz_load(e, res_esyn.raw_flux, sij) for e in measured)
    active = [r for r in rxntext.split(',') if r in res_metab.metabolic_flux and
              r in all_metab_rxns + all_GSM_rxns]
    flux = sum(abs(res_metab.metabolic_flux[r]) for r in active)
    if not measured or not load or flux < min_flux_cutoff: continue
    for rxn in active:
        for enz in measured:
            key, value = rid(rxn, res_metab.metabolic_flux[rxn], enz), mu * flux / load
            kapp[key], kapps_rba[key] = value / 3600, kapp_conversion(value)

# Export kapp summaries and defaults.
exclude = [x for x in open('../exclude_parameterization_list.txt').read().splitlines()
           if x] if os.path.exists('../exclude_parameterization_list.txt') else []
exclude = [enzload_id(x) for x in exclude]
gsm_kapps, gsm_with_data = {}, {}
texts = []
for key, value in kapp.items():
    if key in exclude: continue
    texts.append(f'{key}\t{value}')
    rxn = key.replace('ENZLOAD-', '', 1).rsplit('-', 1)[-2].rsplit('_', 1)[0]
    gsm_kapps.setdefault(rxn, []).append(value)
    if key not in [enzload_id(x) for x in rxns_essential_NP]: gsm_with_data.setdefault(rxn, []).append(value)
write_lines('./kapps_in_vivo.txt', ['rxnid\tkapp (1/s)'] + sorted(texts))
save_csv('./gsm_kapps.csv', gsm_kapps); save_csv('./gsm_kapps_with_prodata.csv', gsm_with_data)
kapp_med = np.median(list(kapp.values())) * 3600 + tol
print('Median kapp (1/h):', kapp_med)
for rxn in rxns_essential_NP + list(enz_rxns_without_enzload): 
    for kapp_key in all_rxn_enzload_pairs.get(rxn, []):
        if kapp_key not in kapps_rba: kapps_rba[kapp_key] = kapp_med; print(f'Using default kapp for {kapp_key} since no kapp could be calculated.')
write_lines('./kapps_per_hr_without_unused_rxns.txt', ['/'] +
            [f"'{k}'\t{v}" for k, v in sorted(kapps_rba.items())] + ['/'])
for rxn in all_metab_rxns:
    tag, base, _, enz = extract_details_from_rxnid(rxn)
    for kapp_key in all_rxn_enzload_pairs.get(rxn, []):
        kapps_rba.setdefault(kapp_key, min(gsm_kapps.get(base, []) + [kapp_med / 3600]) * 3600)
write_lines('./kapps_RBA.txt', ['/'] + [f"'{k}'\t{v}" for k, v in sorted(kapps_rba.items())] + ['/'])
# test if the kapps are feasible in the RBA model
os.system('cd min_flux_violation; module load gams; gams test_kapp.gms --nscale=' + str(nscale2) + output_redirect_str)
try:
    stop_if_run_failed('./min_flux_violation/test_kapp.modelStat.txt')
except Exception:
    # if round_up_default is False, rerun with round_up_default set to True, to see if issue is due to rounding up of kapps
    if not round_up_default and round_up_if_infeas:
        print('RBA run failed. Rerunning with round_up_default set to True to see if issue is due to rounding up of kapps.')
        os.system('python3 C1_calculate_kapp.py 1')
    raise
