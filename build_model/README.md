# Building an RBA Model
- [ ] Build or download the genome-scale metabolic (GSM) model of your chosen organism. Make sure it has been adequately tested beforehand (e.g., via MEMOTE). Put it in the [input](input) folder.
  - **If multiple genes encode a specific protein**: it's recommended to remove all but 1 of them, as these can result in enzyme stoichiometries that lack certain subunits. Prioritizing ones that are less likely to correspond to a target protein (e.g., lower e-value) is ideal.
  - [ ] Make a backup copy of the GSM model if you don't already have one. The model building scripts will add new annotations where possible, to identify potential cofactors.
- [ ] Update [options.yml](input/options.yml) with info requested in comments (i.e., lines starting with `#`).
- [ ] Determine which metabolites the model can't produce, in case any cofactors can't be made. Could use [00_model_reduction.ipynb](00_model_reduction.ipynb) to do this; if the model's default medium isn't accurate to your RBA model's use cases, then update its FVA conditions to include any nutrients that the RBA model may use (by updating the lower bounds on their respective exchange reactions), while excluding all that won't be used in any conditions.
- [ ] Assemble the PROTEIN_stoich_curation.tsv sheet using [01_extract_enzymes_info_uniprot.ipynb](01_extract_enzymes_info_uniprot.ipynb); set `include_unmapped` to False to help with curation of protein compartments, or leave as True to make it faster but less accurate. Follow the instructions at the top of the file. 
  - **How to do this without 01_extract_enzymes_info_uniprot.ipynb** (e.g., if your proteome isn't on UniProt): Using Uniprot or another database, find all genes in your model (and, if alternative splicing is a factor, all proteins that arise from it) that encode enzymes (or subunits of them), transporters, or ribosomes; then, find all subcellular compartments (e.g., mitochondria, ER, peroxisome) where they are present. Create versions of each protein for each compartment they appear in, noting their compartments under subloc_assigned with your model's notation for that compartment (e.g., if "_m" represents the mitochondria, put "m" in subloc_assigned for proteins translated there). Assign their cofactors in the `PTM_stoich_edits` column. Input their names, Uniprot IDs, amino acid sequences, and the compartment they are translated in (under `translation_loc`).
  - **Tips for interpreting/editing the protein stoichiometry:**
    - `trans_stoich_edits`: metabolites to add/subtract from the translation rxn for that protein (e.g., if you want translation to require 2 additional atp_c and release an extra h2o_c, then this section should be "atp_c:2,h2o_c:-1")
    - `PTM_stoich_edits`: like trans_stoich_edits, but for post-translational modifications (PTMs) that happen in a separate rxn. Best spot to put cofactors and PTMs (unless we have evidence that they are added during translation), since it's possible the protein can exist in the cell without them. 
      - "|" characters represent different options for PTMs (e.g., a protein with "fe2_c:1,mg2_c:1|mg2_c:2" as its PTM_stoich_edits value can bind fe2_c and mg2_c once, or mg2_c twice)
      - To let PTMs be applied in compartments other than a protein's translation_loc, can format PTM_stoich_edits entries as a Python dict, with keys representing each compartment (except 'default', which covers all compartments not specified in the dict) and values in the usual PTM_stoich_edits format.
        - Example: For a protein with subloc_assigned set to "c,m" (cytosol, mitochondrion), translation_loc set to "c" (cytosol), and PTM_stoich_edits set to "{'default': 'hemeA_c:1', 'm': 'hemeA_c:1|hemeO_m:1'}", it can bind hemeA_c no matter where it localizes to (since hemeA_c could attach after translation in the cytosol), or bind hemeO_m if it enters the mitochondrion. Useful for representing PTMs that are outside of translation_loc but still usable.
    - Can add [input/01a_organism-specific_changes.py](input/01a_organism-specific_changes.py) to aid in automatic curation after 01_extract_enzymes_info_uniprot.ipynb (e.g., if needing to rerun multiple times while applying changes).
  - [ ] Once [01_extract_enzymes_info_uniprot.ipynb](01_extract_enzymes_info_uniprot.ipynb) is done running, check for any errors or warnings at the bottom.
  - [ ] Check the `comments` column of [PROTEIN_stoich_curation.tsv](input/PROTEIN_stoich_curation.tsv), which may include notes about missing cofactors to consider adding to the model. It's generally recommended to ignore iron-sulfur clusters and (R)-lipoate if they're not found in your GSM model, since finding the right pathways for them can be complex. Some compounds may be listed as cofactors on UniProt despite insufficient evidence that they appear in a given organism (e.g., L-ascorbate in yeasts).
    - **Interpreting comments**: 
      - "Cofactor may be redundant, as it appears in rxns" statements are automatically added to help users avoid potentially overestimating cofactor demand. Review the rxns listed and determine if the compound's usage as a cofactor is accounted for there, as opposed to just being a substrate.
      - "One cluster is coordinated with 3 cysteines and an exchangeable S-adenosyl-L-methionine" statements don't suggest that S-adenosyl-L-methionine is needed, since this statement appears in lipoyl synthase enzymes that cleave S-adenosyl-L-methionine each rxn cycle (https://doi.org/10.1021/bi049528x). If they were other options for cofactors, they would have their own "Name=" statements, instead of just being part of the "Note=" statement for another cofactor.
      - Sometimes, "heme c" is listed as a cofactor despite its iron being included in the "Catalytic activity" subsection (e.g., UniProt accession A0A099P2Y3). In these cases, ferricytochrome c may be a suitable substitute; and if it's accounted for in all of a protein's rxns in the GSM model, then it should be left out of the PTM_stoich_edits column, as with all other redundant cofactors.
- [ ] If your model doesn't follow BiGG naming conventions for tRNAs (e.g., if Arg-tRNA doesn't use ID of "argtrna_c" and its uncharged version's ID isn't "trnaarg_c") or has noncanonical amino acids in its genetic code, then edit PROTEIN_amino_acid_map.txt to use the correct naming schemes. If your model isn't designed for pH 7.2, then also change the 'MW' (molecular weight) column yourself to reflect the actual amino acid MWs.
- Curate the enzymes as you see fit:
  - The newer (recommended) format is stored in [ENZYME_stoich.tsv](input/ENZYME_stoich.tsv).
    - Each row corresponds to a specific enzyme, with a comma-separated list of reactions it catalyzes.
    - By default, the `protein_stoich` column uses gene_src or uniprot IDs from PROTEIN_stoich_curation.tsv. Enzymes will be automatically converted into compartment-specific copies, with their own ENZSYN rxns. Each of these copies will automatically be assembled from the post-translationally modified versions of all proteins in its respective compartment, unless the `protein_stoich` specifies otherwise (e.g., adding a version of a protein from a specific compartment, which may help for transporters).
    - To specify compartments for subunits that don't appear in the same compartment where a rxn takes place (e.g., membrane-bound subunits):
      - If the subunit is used for 1 rxn, then just add the compartment to the `protein_stoich` column (e.g., to specify that protein01 is in compartment "mm", the `protein_stoich` `protein01:1` should be written as `protein01_mm:1`)
      - If it's used for multiple rxns, then it's recommended to use the `compartment-protein_map`, which is a Python dictionary where each key represents a compartment, and is associated with another dictionary containing the protein name and compartment-specific protein ID.
        - Example: To specify that protein01 is in compartment "mm" for rxns in compartment "m", the `compartment-protein_map` column could be set to `{"m":{"protein01":"protein01_mm"}}`. To also add that protein01 should be in "vm" for rxns in "v", then it would become `{"m":{"protein01":"protein01_mm"}, "v":{"protein01":"protein01_vm"}}`.
  - The older format (included for backwards compatibility) is represented by [ENZYME_rxn_pairs.tsv](input/ENZYME_rxn_pairs.tsv); to use it, rename to ENZYME_stoich_curation.tsv.
    - In ENZYME_stoich_curation.tsv, make a row with each reaction's name (`rxn_src` column). If that reaction occurs spontaneously or the proteins needed for it are unknown, put `SPONT` or `UNKNOWN` respectively in its `enz` column and `zeroCost` in its `protein_stoich` column. Otherwise, make a row for every combination of proteins (and how many of them are needed) that can catalyze it (`protein_stoich` column) and assign each row a unique enzyme name (`enz` column).
  - Regardless which format you use, keep in mind that these files are generated using the following assumptions:
      - To reduce the risk of overestimating or inconsistently labeling enzyme demands, when the quantity of a subunit is unknown, assume only 1 copy unless evidence suggests otherwise.
      - By default, all isozymes for a rxn that use a subset of another isozyme's subunits are ignored, since RBA models will generally favor using as few subunits as possible. (e.g., if rxn rxn1_c has the GPR "geneA or (geneA and geneB)", then the GPR will effectively be treated as if it were "geneA")
- **Decide whether to add ribosomes.**
  - **Some pros include:** Avoiding unrealistically high levels of translation, allowing production of ribosome subunits to adjust to what the cell needs instead of remaining at a fixed mass fraction of the proteome.
  - **Some cons include:** Risking inaccurate assignment of ribosome composition or elongation rates; potential impacts can include:
    - Inaccurately treating some ribosome subunits as essential if others are interchangeable with them *in vivo* but not in the model.
  - **If adding them:**
    - Copy [RNA_stoich.xlsx](input/examples/RNA_stoich.xlsx) into `input` folder, fill out a row for each RNA molecule in the model with its name ("RNAid") and sequence. All other info is optional. https://rnacentral.org can be a helpful place to look for RNA info. Remove any placeholder rows.
    - In [input/ribosome_composition.yaml](input/ribosome_composition.yaml), follow the instructions present in the file.
- [ ] To help convert biomass into the RBA model's format, run [01b_biomass-PTM_reconciliation.ipynb](01b_biomass-PTM_reconciliation.ipynb)
  - **If building a new model**:
    - Review biomass composition under [input/biomass_composition.yaml](input/biomass_composition.yaml).
  - **If updating a model that has BIOMASS_RBA.xlsx in its input folder**:
    - Update the cells highlighted in yellow - or if the cell's text is bold, the cells beneath them.
- [ ] Update [phenotype_default.txt](/GAMS/phenotype_default.txt). Add NGAM if not already covered under [SM_rxn_bounds.txt](/GAMS/model/SM_rxn_bounds.txt).
- [ ] Run [02_build_model.ipynb](02_build_model.ipynb)
- [ ] Under [/GAMS/default_paths.txt](/GAMS/default_paths.txt), update the file listed under media_path to match the medium, and rename if needed
- [ ] If you don't have CPLEX access, delete the line `$setGlobal solver cplex` from [default_solver.gms](/GAMS/default_solver.gms).
## $k_{app}$ calculations
- Crucial for defining relationship between rxn flux and enzyme synthesis levels.
- Do this step after determining all other parameters/constraints, to ensure values can be easily reproduced.
- For each set of conditions you're calculating $k_{app}$ values for:
  - [ ] Under [../parameterization/kapp/datasets](../parameterization/kapp/datasets), copy the `(template)` folder and rename however you like. In this new folder:
    - [ ] Adjust `kapp_options.py` to match your organism and paths to your datasets (proteomics - and fluxomics, if available).
    - [ ] run `kapp_calc_pt1.sh`; review sm_j_lumped_mappings.txt to determine if any reaction names differ between your model and flux data, or if any reactions in your flux data are portions of reactions in your model.
      - [ ] Remove the ATP costs of translation from the model's GAM (growth-associated maintenance), if you haven't already.
    - **If you didn't determine ribosome composition beforehand:** 
      - [ ] run `kapp_calc_pt2.sh` and then rerun [02_build_model.ipynb](02_build_model.ipynb) to apply the automatically inferred ribosome composition and kinetics. If changing the ribosome composition from the automatically generated one, run `b2c1.sh` to make sure this doesn't adversely impact its calculations and ensure other ribosome-dependent calculations are performed correctly.
    - **Otherwise:** 
      - [ ] Run `kapp_calc_pt2.sh`; check output for errors. 
    - **In both cases:**
      - If the program stops because the model is infeasible, consider increasing the `tol` value in `min_flux_violation/kapp_calc_defaults.gms`.
      - If the program encounters unspecified errors on SoPlex, or takes too much time/memory, try increasing `iterlim` in [default_solver.gms](/GAMS/default_solver.gms) or using one of the many FVA-based tools provided in `min_flux_violation` (e.g., `FVA_enz_alloc2.gms` in combination with `B2a_kapp_bounds.py`) to estimate and enforce bounds.
    - [ ] Check `min_flux_violation/min_flux_violation.prosynSlack.txt` for any slack variables applied to protein production. These variables are meant to account for cases where a protein's measurements may be inaccurate (e.g., due to excessive rounding of numbers) or paralogs of it may have not been identified. The solver tries to minimize their usage, but checking your dataset for other proteins that should be modeled would also help. Infeasibility may stem from such cases (e.g., ribosome subunits with unidentified paralogs).
    - [ ] Update `kapp_settings.txt`. If additional settings needed, enforce them in `phenotype.txt`.
    - If using these conditions as the default for other calculations (e.g., using these $k_{app}$ values for max yields), update `pheno_root_path` in [phenotype_default_paths.txt](/GAMS/phenotype_default_paths.txt) to the path of this folder.
    - To further improve model quality, it's recommended that you look at `enz_alloc.kappMaxSlack.txt`. If any reactions are listed, it may suggest that some of their enzyme subunits are missing from the proteomics data, and you should consider removing them, looking for homologs in the proteomics data and allowing them to be used instead, or trying other tweaks (e.g., if the rxn is a composite of multiple separate rxns - something you may be able to find out by checking its KEGG page - it may be worth splitting them up into those rxns and reassigning their GPRs accordingly if there's sufficient evidence that not all the subunits are needed to catalyze some of the rxns).
# Decompartmentalization
For simulating relocation of rxns, enzymes, PTMs, and proteins between compartments. 

**If creating a new decompartmentalized model:**
  1. Set `decomp_compartments` value in [input/options.yml](input/options.yml) (e.g., `decomp_compartments: {'m':'c', 'mm':'c'}` allows relocation of rxns/enzymes/PTMs/proteins from compartment `m` to compartment `c`, and from compartment `mm` to compartment `c`).
  2. Make a folder to back up original input files (e.g., `non-decomp` in pkRBA) of `PROTEIN_stoich_curation.tsv`, `ENZYME_stoich.tsv`, and your chosen GSM model.
  3. Run [decompartmentalize.ipynb](decompartmentalize.ipynb). If having problems accessing UniProt, can skip "RBA model updates" heading and run "Rebuild model, update kapps" (though this won't update cofactor decompartmentalization as thoroughly as a full run).
  4. Make a backup of the same input files you backed up in Step 2 (e.g., putting them in `decompartmentalized` folder).

**If replicating results for a model that has already been decompartmentalized:**
  1. Copy its backup input files (e.g., in `decompartmentalized` folder) into [input/](input/).

Then, rerun [02_build_model.ipynb](02_build_model.ipynb).
## Reversing model decompartmentalization
- **Best method for keeping up-to-date models**: run [decompartmentalization_reversal.ipynb](decompartmentalization_reversal.ipynb).
- **Best method for replicating prior results**: Copy the non-decompartmentalized input files (e.g., in `non-decomp` folder) into [input/](input/).
  - However, database updates made since the non-decompartmentalized input files were generated won't be applied (e.g., if UniProt enzyme stoichiometry predictions were updated when creating the decompartmentalized model, they won't be included in the non-decompartmentalized version).