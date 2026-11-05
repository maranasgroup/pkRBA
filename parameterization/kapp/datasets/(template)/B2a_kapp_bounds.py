# Should be run after min_flux_violation.gms (but most useful when run before enz_alloc.gms, so consider separating these into two separate scripts or wedging this between them in B2)
# help find maximum possible kapp, using FVA_enz_alloc2.txt (table with 1+ spaces separating values) to reference the min and max fluxes for all rxns
from kapp_options import *
from collections import defaultdict
EMPTY = {"/", "\n", "/\n"}

def read_lines(path):
    with open(path) as f:
        return [line for line in f if line.strip() not in EMPTY]

# Reaction/enzyme lookup tables.
rxn_enzload_map = defaultdict(list)
enz_load_map = defaultdict(list)
enzsyn_load_map = defaultdict(list)
enzload_syn_map = defaultdict(list)
proteins_and_locations = {}
prosyn_and_proteins = {}

for line in read_lines(os.path.join(path_model, "proteins_and_locations.txt")):
    protein, location = line.split("'.(")
    protein = protein.lstrip("'")
    locations = [
        item.strip().strip("'\"")
        for item in location.strip().strip("()").split(",")
        if item.strip()
    ]
    proteins_and_locations[protein] = locations
    for rxn in locations:
        prosyn_and_proteins[rxn] = protein

v_exp_pro = {}
for line in read_lines("v_exp_pro.txt"):
    met_id, value = line.replace("'", "").split()
    v_exp_pro[met_id] = float(value)

mets, s = set(), {}
met_to_producers = defaultdict(list)
mets_freely_exchangeable = {}

with open(sij_path) as f:
    for line in f:
        if not line.strip() or line.strip() in EMPTY:
            continue

        met_id, rxn_id, coefficient = (
            line.replace("'.'", "' '").replace("'", "").split()
        )
        coefficient = float(coefficient)
        mets.add(met_id)
        s.setdefault(rxn_id, {})[met_id] = coefficient

        if coefficient > 0:
            met_to_producers[met_id].append(rxn_id)

        if rxn_id.startswith("ENZLOAD-"):
            _, base, direction, _ = extract_details_from_rxnid(rxn_id)
            reaction_id = f"RXN-{base}_{direction}"
            rxn_enzload_map[reaction_id].append(rxn_id)
            enz_load_map[met_id].append(rxn_id)

def read_bounds(path, extra_columns=[]):
    with open(path) as file:
        rows = file.read().splitlines()
    header = extra_columns + rows[0].split()
    return pd.DataFrame([row.split() for row in rows[1:] if row], columns=header).set_index("rxn")

df_bounds_calc = read_bounds(os.path.join("min_flux_violation", "FVA_enz_alloc2.txt"))
df_bounds = read_bounds(os.path.join("min_flux_violation", "FVA_enz_alloc2.bounds_all.txt"), ["rxn"])

for column in ("min_flux", "max_flux", "low_stat", "high_stat"):
    df_bounds_calc[column] = pd.to_numeric(
        df_bounds_calc[column], errors="coerce"
    )

for column in ("lo", "up"):
    df_bounds[column] = pd.to_numeric(df_bounds[column], errors="coerce")

df_bounds[["lo", "up"]] /= nscale2

# Apply model bounds where statistical bounds are unavailable.
common = df_bounds_calc.index.intersection(df_bounds.index)
mask = df_bounds_calc.loc[common, "low_stat"] != 1
df_bounds_calc.loc[common[mask], "min_flux"] = df_bounds.loc[common[mask], "lo"]

mask = df_bounds_calc.loc[common, "high_stat"] != 1
df_bounds_calc.loc[common[mask], "max_flux"] = df_bounds.loc[common[mask], "up"]

# Add reactions missing from the calculated bounds.
missing_rxns = s.keys() - df_bounds_calc.index

if missing_rxns:
    missing = pd.DataFrame(index=list(missing_rxns))
    missing["min_flux"] = df_bounds["lo"].reindex(missing.index).fillna(0)
    missing["max_flux"] = df_bounds["up"].reindex(missing.index).fillna(0)
    missing[["low_stat", "high_stat"]] = 1

    for column in df_bounds_calc.columns:
        if column not in missing:
            missing[column] = np.nan

    df_bounds_calc = pd.concat([df_bounds_calc, missing])

invalid = df_bounds_calc["min_flux"] > df_bounds_calc["max_flux"]

if invalid.any():
    invalid_rxns = df_bounds_calc.index[invalid].tolist()
    df_bounds_calc.loc[invalid, "max_flux"] = df_bounds_calc.loc[invalid, "min_flux"]
    print(f"Warning: min_flux exceeded max_flux; set max_flux equal to min_flux for {invalid_rxns}")

df_bounds_calc["min_flux_max_enz"] = df_bounds_calc["min_flux"]
df_bounds_calc["max_flux_max_enz"] = df_bounds_calc["max_flux"]

def max_enz_flux(rxn):
    return df_bounds_calc["max_flux_max_enz"].get(rxn, 0)

product_consumers = defaultdict(list)

for reaction, stoichiometry in s.items():
    for metabolite, coefficient in stoichiometry.items():
        if coefficient < 0:
            product_consumers[metabolite].append(reaction)
    if reaction.startswith('RXN-') and len(stoichiometry) == 2 and df_bounds_calc.at[reaction, "max_flux_max_enz"] > 0:
        _, base, _, _ = extract_details_from_rxnid(reaction)
        # if "SBO:0000176" in annotations in model, then ignore it since it's not a transport reaction
        if model.reactions.get_by_id(base).__dict__.get("_annotation", {}).get("sbo") == "SBO:0000176":
            continue
        # if the max coefficient is the inverse of the min coefficient, then it's a transport reaction and we can add the positive metabolite to the list of freely exchangeable metabolites
        if max(stoichiometry.values()) == -min(stoichiometry.values()):
            # add the positive metabolite to the list of freely exchangeable metabolites alongside its negative metabolite (key) and the reaction (value) that connects them
            positive_metabolite = next(met for met, coeff in stoichiometry.items() if coeff > 0)
            negative_metabolite = next(met for met, coeff in stoichiometry.items() if coeff < 0)
            mets_freely_exchangeable.setdefault(negative_metabolite, {}).setdefault(positive_metabolite, []).append(reaction)

# find all PTM rxns usable for each protein using proteins_and_locations and df_bounds_calc
ptm_rxns_by_protein = defaultdict(list)
proteins_by_ptm_rxn = defaultdict(list)
for protein, locations in proteins_and_locations.items():
    for location in locations:
        key_products = [
            metabolite
            for metabolite, coefficient in s[location].items()
            if coefficient > 0 and metabolite.startswith("PRO-")
        ]
        for product in key_products:
            for consumer in product_consumers.get(product, ()):
                if consumer.startswith("PTM-") and max_enz_flux(consumer) > 0:
                    ptm_rxns_by_protein[protein].append(consumer)
                    proteins_by_ptm_rxn[consumer].append(protein)

# show the metabolites that are freely exchangeable
# for negative_metabolite, positive_metabolites in mets_freely_exchangeable.items():
#     for positive_metabolite, reactions in positive_metabolites.items():
#         print(f"{negative_metabolite} is freely exchangeable with {positive_metabolite} via reactions: {', '.join(reactions)}")

for enzload_rxns in rxn_enzload_map.values():
    for enzload_rxn in enzload_rxns:
        enzyme = next(metabolite for metabolite, coefficient in s[enzload_rxn].items()
            if coefficient < 0
        )
        for enzsyn_rxn in met_to_producers.get(enzyme, ()):
            if enzsyn_rxn.startswith("ENZSYN-"):
                enzload_syn_map[enzload_rxn].append(enzsyn_rxn)
                enzsyn_load_map[enzsyn_rxn].append(enzload_rxn)

proteins_and_locations_used = {}
rxns_with_min_metabolite_allocations = set()

while True:
    changes, allocations = 0, {}

    for rxn in df_bounds_calc.index:
        if rxn.startswith("PROSYN-"):
            waste_rxn = rxn.replace("PROSYN-", "PROWASTE-", 1)
            product_prefix = "PRO-"
        elif rxn.startswith("PTM-"):
            waste_rxn = rxn.replace("PTM-", "PTMWASTE-", 1)
            product_prefix = "PTMPRO-"
        else:
            continue

        key_products = [
            metabolite
            for metabolite, coefficient in s[rxn].items()
            if coefficient > 0 and metabolite.startswith(product_prefix)
        ]
        syn_flux = df_bounds_calc.at[rxn, "min_flux_max_enz"]
        waste_flux = df_bounds_calc["min_flux_max_enz"].get(waste_rxn, 0)

        for product in key_products:
            consumers = [
                consumer
                for consumer in product_consumers.get(product, ())
                if max_enz_flux(consumer) > 0
                and not consumer.startswith(("PROWASTE-", "PTMWASTE-"))
            ]
            if len(consumers) > 1:
                continue
            if len(consumers) == 0:
                if not (rxn.startswith("PTM-") and df_bounds_calc.at[rxn, "min_flux_max_enz"] == 0):
                    continue
                # check if its MET- reactants and products are freely exchangeable w/ metabolites in other usable PTMs for the same protein, and if so, ban this PTM and all rxns using its product
                met_stoich = {metabolite: coefficient for metabolite, coefficient in s[rxn].items() if metabolite.startswith("MET-")}
                proteins = proteins_by_ptm_rxn.get(rxn, [None])
                # print(f"Checking if {rxn} has freely exchangeable metabolites with other PTMs for proteins {proteins}")
                for protein in proteins:
                    for other_ptm in ptm_rxns_by_protein.get(protein, []):
                        if other_ptm == rxn or df_bounds_calc.at[other_ptm, "max_flux_max_enz"] == 0:
                            continue
                        # check if all its metabolites are freely exchangeable w/ this PTM's metabolites
                        other_met_stoich = {metabolite: coefficient for metabolite, coefficient in s[other_ptm].items() if metabolite.startswith("MET-")}
                        match = True if met_stoich == other_met_stoich else False
                        # if not match and all(met in mets_freely_exchangeable and any(other_met in mets_freely_exchangeable[met] for other_met in other_met_reactants + other_met_products) for met in met_reactants + met_products):
                        #     # ban this PTM and all rxns using its product
                        #     print(f"Banning {rxn} and all rxns using its product {product} because it has no consumers and is freely exchangeable with {other_ptm} for protein {protein}")
                        #     df_bounds_calc.loc[[rxn] + product_consumers.get(product, []), ["min_flux_max_enz", "max_flux_max_enz"]] = 0
                        #     changes += 1
                        if match:
                            df_bounds_calc.loc[[rxn] + product_consumers.get(product, []), ["min_flux_max_enz", "max_flux_max_enz"]] = 0
                            changes += 1
                            # remove rxn from proteins_by_ptm_rxn to avoid checking it again
                            proteins_by_ptm_rxn.pop(rxn, None)
                            # print(f"Banning {rxn} and all rxns using its product {product} because it has no consumers and is identical to {other_ptm} for protein {protein}")
                continue

            consumer = consumers[0]
            allocations.setdefault(consumer, {})[product] = max(
                0, (syn_flux - waste_flux) / -s[consumer][product]
            )
            protein = prosyn_and_proteins.get(rxn)
            if protein is not None:
                proteins_and_locations_used.setdefault(protein, set()).add(rxn)

    for protein, used_locations in proteins_and_locations_used.items():
        if len(used_locations) == 1 and protein in v_exp_pro:
            location = next(iter(used_locations))
            new_val = min(
                df_bounds_calc.at[location, "max_flux_max_enz"],
                v_exp_pro[protein],
            )
            if new_val > df_bounds_calc.at[location, "min_flux_max_enz"]:
                df_bounds_calc.at[location, "min_flux_max_enz"] = new_val
                changes += 1

        unused = set(proteins_and_locations[protein]) - used_locations
        for location in unused:
            if df_bounds_calc.at[location, "max_flux_max_enz"] == 0:
                continue

            df_bounds_calc.loc[
                location, ["min_flux_max_enz", "max_flux_max_enz"]
            ] = 0
            changes += 1

            products = [
                metabolite
                for metabolite, coefficient in s[location].items()
                if coefficient > 0 and metabolite.startswith("PRO-")
            ]
            consumers = {
                consumer
                for product in products
                for consumer in product_consumers.get(product, ())
            }
            existing = consumers & set(df_bounds_calc.index)

            if existing:
                df_bounds_calc.loc[
                    list(existing),
                    ["min_flux_max_enz", "max_flux_max_enz"],
                ] = 0
                changes += 1

    for reaction, product_allocations in allocations.items():
        new_value = min(
            max(
                df_bounds_calc.at[reaction, "min_flux_max_enz"],
                min(product_allocations.values()),
            ),
            df_bounds_calc.at[reaction, "max_flux_max_enz"],
        )

        if new_value > df_bounds_calc.at[reaction, "min_flux_max_enz"]:
            df_bounds_calc.at[
                reaction, "min_flux_max_enz"
            ] = new_value
            changes += 1

        enzloads = [enzload for enzload in enzsyn_load_map.get(reaction, ()) if max_enz_flux(enzload) > 0]

        for enzload in enzsyn_load_map.get(reaction, ()):
            value = (
                df_bounds_calc.at[reaction, "min_flux_max_enz"]
                / len(enzloads)
                if enzloads
                else 0
            )
            value = min(value, max_enz_flux(enzload))

            if value > df_bounds_calc.at[enzload, "min_flux_max_enz"]:
                df_bounds_calc.at[
                    enzload, "min_flux_max_enz"
                ] = value
                changes += 1

            _, base, direction, _ = extract_details_from_rxnid(enzload)
            rxns_with_min_metabolite_allocations.add(
                f"RXN-{base}_{direction}"
            )

    if not changes:
        break

for reaction, enzloads in rxn_enzload_map.items():
    max_enzload = sum(max_enz_flux(rxn) for rxn in enzloads)
    min_enzload = sum(
        df_bounds_calc.at[rxn, "min_flux_max_enz"]
        for rxn in enzloads
        if rxn in df_bounds_calc.index
    )
    min_flux = df_bounds_calc["min_flux_max_enz"].get(reaction, 0)
    max_flux = df_bounds_calc["max_flux_max_enz"].get(reaction, 0)

    if min_enzload == 0 and (max_flux > 0 or min_flux > 0):
        print(
            f"Warning: min. ENZLOAD of 0 for {reaction} w/ bounds {(min_flux,max_flux)}. Enzymes: {enzloads}. Has min_metabolite_allocations: {reaction in rxns_with_min_metabolite_allocations}"
        )
    if max_enzload == 0 and min_flux > 0:
        print(
            f"Warning: max. ENZLOAD of 0 for {reaction} w/ bounds {(min_flux,max_flux)}. Enzymes: {enzloads}. Has min_metabolite_allocations: {reaction in rxns_with_min_metabolite_allocations}"
        )
        # restore the enzload fluxes to their original values, since the reaction is still active and has a non-zero min_flux
        for enzload in enzloads:
            if enzload in df_bounds_calc.index:
                df_bounds_calc.at[enzload, "max_flux_max_enz"] = RBA_defaults.vmax
                #TODO: Also find ways to fix other upstream reactions that may have been set to 0 due to this issue, such as ENZSYN and PROSYN reactions

# sort by reaction ID for easier comparison with other files
df_bounds_calc = df_bounds_calc.sort_index()

df_bounds_calc.to_csv(
    "min_flux_violation/enz_alloc.bounds.tsv", sep="\t"
)

with open("min_flux_violation/enz_alloc.bounds.txt", "w") as file:
    max_lengths = {
        column: max(len(str(value)) for value in df_bounds_calc[column].tolist() + [column])
        for column in df_bounds_calc.columns
    }
    max_lengths["rxn"] = max(len(rxn) for rxn in df_bounds_calc.index.tolist() + ["rxn"])

    file.write(
        f"{'   '.ljust(max_lengths['rxn'])} "
        f"{'min_flux_max_enz'.ljust(max_lengths['min_flux_max_enz'])} "
        f"{'max_flux_max_enz'.ljust(max_lengths['max_flux_max_enz'])} "
        f"{'low_stat'.ljust(max_lengths['low_stat'])} "
        f"{'high_stat'.ljust(max_lengths['high_stat'])}\n"
    )

    for rxn, row in df_bounds_calc.iterrows():
        file.write(
            f"{rxn:<{max_lengths['rxn']}} "
            f"{str(row['min_flux_max_enz']).ljust(max_lengths['min_flux_max_enz'])} "
            f"{str(row['max_flux_max_enz']).ljust(max_lengths['max_flux_max_enz'])} "
            f"{str(row['low_stat']).ljust(max_lengths['low_stat'])} "
            f"{str(row['high_stat']).ljust(max_lengths['high_stat'])}\n"
        )

invalid = df_bounds_calc[
    df_bounds_calc["min_flux_max_enz"] > df_bounds_calc["max_flux_max_enz"]
]

if not invalid.empty:
    raise ValueError("Some reactions have min_flux_max_enz greater than max_flux_max_enz:\n"
        f"{invalid.index.tolist()}"
    )