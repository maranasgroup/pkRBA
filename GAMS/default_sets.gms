sets
i
$include %species_path%
pro
$include %pro_path%
AA
$include %aa_path%
tRNA(AA,i)
$include %tRNA_path%
i_enz(i)
j
$include %rxns_path%
g
$include %genes_path%
c
$include %compartments_path%
c_mito(c)
$include %compartments_mito_path%
j_c(j,c)
$include %rxn_compartment_path%
expression(g,j)
$include %expression_path%
sm_j "stoichiometric model (SM) reactions"
$include %sm_rxns_path%
prosyn(j) 
$include %prosyn_path%
ptm(j)
$include %gms_path%RBA_rxns_PTM.txt
pro_prosyn(pro,j)
$include %proteins_and_locations_path%
pro_syn_waste(j,j)
$include %gms_path%pro_syn_waste.txt
ptmwaste(j) "sinks for post-translationally modified proteins"
$include %gms_path%RBA_rxns_ptmwaste.txt
r(*)
$include %ribosomes_path%
r_active(r) "active ribosomes"
pro_ribo(r,pro)
$include %gms_path%pro_ribo_subunits.txt
j_ribo(r,j)
$include %j_ribo_path%
trans(r,j)
$include %translation_path%
trans_by_any_ribo(j)
$include %trans_by_any_ribo_path%
prowaste(j)
$include %prowaste_path%
biowaste(j)
$include %biowaste_path%
uptake(j) "uptake rxns"
$include %uptake_path%
secretion(j) "secretion rxns"
$include %secretion_path%
media(j) "allowed uptake rxns based on simulated media conditions"
$include %media_path%
rxns_biomass(j)
$include %biomass_path%
j_bio_prot(j)
$include %j_bio_prot_path%
enzsyn(j)
$include %enzsyn_path%
enzload(j)
$include %enzload_path%
enzload_rxn_coupling(j,j)
$include "%enzload_rxn_coupling_path%"
rxns_with_enz(j)
$include %rxns_with_enz_path%
;
