* Run FVA (Flux Variability Analysis) on the RBA model
* To use in conjunction with another script, set the objective function to its optimal value from that script and then include this file
$if not setGlobal gms $setGlobal gms %system.FN%

$INLINECOM /*  */
$include ../paths.txt
$offInline

* run RBA without bisection, using default settings and outputs
$if not setGlobal gms $setGlobal gms %system.FN%
$include kapp_calc_defaults.gms

prosynSlackUB.up(pro_measured) = inf;
prosynSlackLB.up(pro_measured) = 1 + %bi_tol%;

Sets
rxns_tested(j) "rxns to include in FVA" /
#prosyn
/
dummy(j) "dummy set for iterating over rxns"
;

Parameters
high(j), low(j)
;

Variables
z_fva
;

Equations
Obj_FVA
;

Obj_FVA..           z_fva =e= sum(j$dummy(j),v(j));

$onMultiR
model rba /all-EnzCap-minPro-ProwasteLim-prot_mass_fraction_max_eqn-membrane_SA_vivo-Weight_Prowaste-kappEstSlack-Obj5-membrane_SA_eqn-v_trans_UB_fx-v_trans_LB_fx-v_trans_AA_sum_def-v_sum_def/;
rba.optfile = 1;
$offMulti

*Updates output while solving, so errors won't get rid of output
file f_fva /%system.FN%.txt/;
put f_fva;
* Clear old file
put "rxn",system.tab,"upper-bound"/;
* Append old file
f_fva.ap=1;

$if exist %phenotype_path% $include %phenotype_path%

alias(j_fva,j);

loop(j_fva$(rxns_tested(j_fva)),
    dummy(j)$(ord(j) ne ord(j_fva)) = no;
    dummy(j)$(ord(j) eq ord(j_fva)) = yes;
    display j_fva;
    
    Solve rba using lp maximizing z_fva;
    If(rba.modelstat eq 1,
        high(j_fva) = z_fva.l;
		low(j_fva) = 0;
		put f_fva;
		put j_fva.tl:0:50,system.tab,high(j_fva):0:8/;
		putclose f_fva
        );
    );

solve rba using lp minimizing z_fva;
$include %model_root_path%GAMS/default_outputs.gms
