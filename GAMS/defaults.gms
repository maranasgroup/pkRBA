* Default setup for RBA runs, including options, sets/parameters, variables, and output file names
* Reduces redundancy in the code and makes sweeping changes easier to implement
* default paths or settings specific to your model should go in phenotype_default_paths.txt or phenotype_default.txt respectively
* usage: $include <path_to_this_file>

$INLINECOM /*  */
* use default settings if not set already
$if not setGlobal gms $setGlobal gms %system.FN% 
$if not setGlobal model_root_path $setGlobal model_root_path ../ 
$if not setGlobal gms_path $setGlobal gms_path %system.FP%model/
$if not setGlobal settings_path $setGlobal settings_path %system.FP%default_paths.txt
$include %settings_path%
* default model-/phenotype-specific settings; not needed, so only loaded if it exists
$if exist %phenotype_default_paths_path% $include %phenotype_default_paths_path%
* Scale values of all variables by a factor, then when write to file, descale them
$setGlobal vmax 1e3 /* upper bound for fluxes */
$setGlobal vmin 1e-12 /* lowest positive value distinguishable from 0 by solver */
$setGlobal nscale 1e5 /* multiplier for fluxes, to reduce solver precision issues */
$setGlobal mu_exp 0.1 /* default measured growth rate used as placeholder (e.g., for yields) */
$setGlobal mu %mu_exp% /* default growth rate used as placeholder (e.g., for yields) */

* for performing bisection
** determines tolerance (i.e., stops when upper/lower bounds are negligibly far apart)
$setGlobal bi_tol 1e-4
** max number of times bisection occurs
$setGlobal max_iter 100

* Default: ≤100% of protein mass can be in mitochondria or its membrane
$setGlobal max_allowed_mito_proteome_allo_fraction 1
* Enforce part of proteome allocate to non-modeled protein
$setGlobal nonmodeled_proteome_allocation 0

* default solver settings
$include %system.FP%default_solver.gms

* remove existing modelStat file to avoid issues w/ model status detection       
file ff /%gms%.modelStat.txt/; putclose ff '';

$include %system.FP%default_sets.gms
r_active(r)$(sum(pro$pro_ribo(r,pro), 1) gt 0) = yes;
alias (j1,j);

Parameters
S(i,j)
$include %sij_path%
dir(sm_j,j) "SM rxn, RBA rxn, and direction (-1 if RBA rxn is the reverse of SM rxn, 1 otherwise)"
$include "%sm_rxn_pairs_path%"
v_trans_UB(AA) "max molar fraction of all AAs that each AA can comprise"
$ifThenE setGlobal v_trans_UB_AA_path
$include %v_trans_UB_AA_path%
$else
/#AA 1/
$endIf
v_trans_LB(AA) "same as above, but min instead" /#AA 0/
NAA(j)
$include %prolen_path%
mv_exp_prot "Total flux thru BIOSYN-PROTTOBIO (i.e., total flux of rxn incorporating protein mass into biomass)"
$if exist %prot_mass_flux_path% $include %prot_mass_flux_path%
kapp(j)
$if not %ignore_kapps%==1 $include %kapp_path%
kribo(r)
$include %kribo_path%
bi_UB "upper bound for bisection" /1/
bi_LB "lower bound for bisection" /0/
max_iter /%max_iter%/
mu
prowaste_predicted_min(j)/
$ifThenE card(pro_ribo) == 0
$if exist %prowaste_max_from_pro_with_data_path% $include %prowaste_max_from_pro_with_data_path%
$if exist %prowaste_max_from_pro_without_data_path% $include %prowaste_max_from_pro_without_data_path%
$else
$if exist %prowaste_min_from_pro_with_data_path% $include %prowaste_min_from_pro_with_data_path%
$if exist %prowaste_min_from_pro_without_data_path% $include %prowaste_min_from_pro_without_data_path%
$endIf
/
prowaste_predicted_max(j)/
$if exist %prowaste_max_from_pro_with_data_path% $include %prowaste_max_from_pro_with_data_path%
$if exist %prowaste_max_from_pro_without_data_path% $include %prowaste_max_from_pro_without_data_path%
/
prot_mass_fraction_max(j) 
$if exist %prot_mass_fraction_max_path% $include %prot_mass_fraction_max_path%
use_prot_mass_fraction_max_eqn(j) "1 if prot_mass_fraction_max_eqn(j) should be used, 0 otherwise" / #j 1 /
use_membrane_SA_vivo(c) "1 if membrane_SA_vivo(c) equations should be used, 0 otherwise" / #c 0 /
membrane_SA_exp(c) "in vivo membrane surface area for each compartment, for membrane protein allocation"
$if exist %membrane_SA_exp_path% $include %membrane_SA_exp_path%
beta(c) "growth rate dependent change in surface area of compartment's membrane due to change in SA:V"
fMSA(c) "fraction of membrane surface area available for proteins" / #c 1 /
delta(c) "nonideality in surface area usage including misfolded proteins or misassembled protein complexes" / #c 0 /
c1(c) "slope of linear relationship between compartment SA:V and growth rate, set to lowest E. coli value by default" / #c -1.384 /
theta(c) "intercept of linear relationship between compartment SA:V and growth rate, set to lowest E. coli value by default" / #c 0.048 /
SA(j) "membrane surface area per unit flux for each reaction, for membrane protein allocation"
cap_gross_flux(sm_j) "whether gross flux is limited by net flux bounds for a given rxn; no by default" / #sm_j 0 /
minPro_fraction / 1 /
;
$include %protein_membrane_surface_area_path%
beta(c) = (fMSA(c) - delta(c)) * (c1(c) / theta(c));

i_enz(i) = yes$(sum(j$(enzsyn(j) and S(i,j) gt 0),1) gt 0);

Variables
z, v(j), v_sum, v_trans(AA), v_trans_AA_sum, membrane_SA(c)
z_dummy "dummy objective used just to test for feasibility";
Table v_sm(sm_j, *)
$include %sm_rxn_bounds_path%
;
Parameters
v_sm_lo(sm_j) "lower bounds for fluxes represented by sm_j"
v_sm_up(sm_j) "upper bounds for fluxes represented by sm_j"
;
v_sm_up(sm_j) = v_sm(sm_j, 'v_sm_up');
v_sm_lo(sm_j) = v_sm(sm_j, 'v_sm_lo');

*** FLUX LOWER AND UPPER BOUNDS (unless otherwise needed, use v_sm_lo and v_sm_up to alter aggregate flux bounds) ***
v.lo(j) = 0; v.up(j) = %vmax%; 

* Disable all biomass reactions, to avoid accidental usage
* User should activate condition-specific biomass reaction in phenotype.txt file
v.fx(j)$rxns_biomass(j) = 0;
* Disable all biowaste reactions (meant for 0-growth simulations), to avoid accidental usage
v.fx(j)$biowaste(j) = 0;

* optional model-/phenotype-specific settings 
$if exist %phenotype_default_path% $include %phenotype_default_path%
* optional task-specific phenotype file
$if exist %phenotype_path% $include %phenotype_path%
*allow bisection
mu=%mu%;

* scale all bounds by nscale, to improve solver precision; any v bounds defined after here should follow suit
v_sm_up(sm_j) = v_sm_up(sm_j) * %nscale%;
v_sm_lo(sm_j) = v_sm_lo(sm_j) * %nscale%;
v.up(j) = v.up(j) * %nscale%;
v.lo(j) = v.lo(j) * %nscale%;

Equations
Obj, z_dummy_def, Stoic, v_sum_def, v_trans_def, v_trans_AA_sum_def, v_trans_UB_fx, v_trans_LB_fx, SM_bounds_LB, SM_bounds_UB, SM_gross_flux, RiboCapacity, RiboCapacity2, FlexibleRiboCapacity, minPro, ProwasteLim, prot_mass_fraction_max_eqn, membrane_SA_eqn, membrane_SA_vivo(c), NonModelProtAllo, MitoProtAllo, EnzCap(j)
;

* Should be able to reflect MitoProtAllo and separate mm intermembrane proteins from mm transmembrane proteins
* ProtAllo(c)$(max_prot_mass_fraction(c)).. v('BIOSYN-PROT_'c) - v('BIOSYN-PROT_IS_'c) =l= v('BIOSYN-PROTMODELED') * max_prot_mass_fraction(c);
* MitoProtAllo.. v('BIOSYN-PROT_m') + v('BIOSYN-PROT_mm') =l= %max_allowed_mito_proteome_allo_fraction% * v('BIOSYN-PROTMODELED');

v_sum_def..		v_sum =e= sum(j, v(j));
Obj..			z =e= v_sum;
z_dummy_def..	z_dummy =e= 0;
Stoic(i)..		sum(j, S(i,j)*v(j)) =e= 0;
v_trans_def(AA).. v_trans(AA) =e= sum(i$tRNA(AA,i),sum(j$(S(i,j) lt 0 and (sum(r,trans(r,j)) or trans_by_any_ribo(j))),v(j)*(-S(i,j))));
v_trans_AA_sum_def..	v_trans_AA_sum =e= sum(AA,v_trans(AA));
v_trans_UB_fx(AA).. v_trans(AA) =l= v_trans_AA_sum * v_trans_UB(AA);
v_trans_LB_fx(AA).. v_trans(AA) =g= v_trans_AA_sum * v_trans_LB(AA);
SM_bounds_LB(sm_j)..	v_sm_lo(sm_j) =l= sum(j$dir(sm_j,j), dir(sm_j,j)*v(j));
SM_bounds_UB(sm_j)..	sum(j$dir(sm_j,j), dir(sm_j,j)*v(j)) =l= v_sm_up(sm_j);
SM_gross_flux(sm_j).. cap_gross_flux(sm_j)*sum(j$dir(sm_j,j), v(j)) =l= max(abs(v_sm_up(sm_j)), abs(v_sm_lo(sm_j)));
RiboCapacity(r_active(r))$(card(trans_by_any_ribo) eq 0).. sum(j$j_ribo(r,j),v(j)) =e= mu / kribo(r) * sum(j$trans(r,j), NAA(j) * v(j));
* Alternate version added when there are any trans_by_any_ribo, to allow for more flexibility in ribosome usage (prior version used to reduce solution degeneracy)
RiboCapacity2(r_active(r))$(card(trans_by_any_ribo) > 0).. sum(j$j_ribo(r,j),v(j)) =g= mu / kribo(r) * sum(j$trans(r,j), NAA(j) * v(j));
* FlexibleRiboCapacity only included if there are any trans_by_any_ribo
FlexibleRiboCapacity$(card(trans_by_any_ribo) > 0).. sum(r_active(r),sum(j$j_ribo(r,j),v(j) * kribo(r))) =g= mu * (sum(r_active(r),sum(j$trans(r,j), NAA(j) * v(j))) + sum(j$trans_by_any_ribo(j), NAA(j) * v(j)));

* Also known as "minimum protein from sink (MPFS)" constraints, these require a protein's mass fraction of the proteome to at least match its predicted value from the data. Exceptions are made for proteins that the user knocks out, to allow testing of these knockouts.
minPro(j)$(v.up(j) ne 0 and prowaste_predicted_min(j) ne 0).. v(j) * mv_exp_prot =g= v('BIOSYN-PROTTOBIO') * prowaste_predicted_min(j) * minPro_fraction;

* Limits PROWASTE to the predicted value (if available); to allow any waste value, make the LHS >= 0. Values kept relative to total protein mass from kapp calculations for reliable scaling across different conditions.
ProwasteLim(j)$(prosyn(j) or ptm(j)).. sum(j1$pro_syn_waste(j,j1),v(j1)) * mv_exp_prot =l= v('BIOSYN-PROTTOBIO') * prowaste_predicted_max(j);
* Compartmental proteome mass (CPM) constraints. Sets upper limits on the mass fraction of the proteome that each protein type (e.g., membrane-bound cytoplasmic, soluble mitochondrial, etc.) can comprise, based on the predicted proteome allocation from the data
prot_mass_fraction_max_eqn(j)$(prot_mass_fraction_max(j)).. use_prot_mass_fraction_max_eqn(j)*v(j) =l= (prot_mass_fraction_max(j) / mv_exp_prot) * v('BIOSYN-PROTTOBIO');

* NonModelProtAllo and MitoProtAllo are included only for backwards compatibility with scRBA and rtRBA
NonModelProtAllo(j)$sameAs(j,'BIOSYN-PROTMODELED')..	v(j) =e= (1 - %nonmodeled_proteome_allocation%) * sum(j1$sameas(j1,'BIOSYN-PROTTOBIO'),v(j1));
MitoProtAllo(j)$sameAs(j,'BIOSYN-PROTMODELED')..		sum(c_mito, sum(j_bio_prot$(j_c(j_bio_prot,c_mito)),v(j_bio_prot))) =l= %max_allowed_mito_proteome_allo_fraction% * v(j);

membrane_SA_eqn(c)..	membrane_SA(c) =e= sum(j$(SA(j) and j_c(j,c)), SA(j) * v(j));
membrane_SA_vivo(c)$membrane_SA_exp(c)..	use_membrane_SA_vivo(c)*(-beta(c)*mu + (membrane_SA(c)/%nscale%)) =l= use_membrane_SA_vivo(c)*(-beta(c)*%mu_exp% + membrane_SA_exp(c));
EnzCap(j)$rxns_with_enz(j).. v(j) * mu =e= sum(j1$enzload_rxn_coupling(j1,j),v(j1) * kapp(j1))

*** BUILD OPTIMIZATION MODELS ***
Model 
rba /all/
fba /rba-RiboCapacity-RiboCapacity2-FlexibleRiboCapacity-EnzCap/
;
rba.optfile = 1;
fba.optfile = 1;

file fi_v /%gms%.flux.txt/;
file fi_vgs /%gms%.flux_gamsscaled.txt/;
