* Run FVA (Flux Variability Analysis) on the RBA model
* To use in conjunction with another script, set the objective function to its optimal value from that script and then include this file
$if not setGlobal gms $setGlobal gms %system.FN%

$setGlobal ignore_blocked_j 1
$setGlobal rerun 0

$INLINECOM /*  */
$include ../paths.txt
$offInline

* run RBA without bisection, using default settings and outputs
$include kapp_calc_defaults.gms

prosynSlackUB.up(pro_measured) = inf;
prosynSlackLB.up(pro_measured) = 1 + %bi_tol%;

$onMultiR
Sets
rxns_to_test(j)
$ifThen %rerun%==1
$include "rxns_to_test.txt"
$else $include %rxns_path%
$endIf
blocked(j) //
;
scalar nonzero_ub
$offMulti

loop(enzload_rxn_coupling(j,j1),
    if (v.up(j1) eq 0,
        v.up(j) = 0;
    );
);

$include %model_root_path%find_blocked_rxns.gms
blocked(j)$(v.lo(j) gt -%vmin% and v.up(j) lt %vmin%)=yes;

variables z_fva;

equation obj_fva; obj_fva.. z_fva =e= sum(j$rxns_to_test(j),v(j));
$onMultiR
model rba /all-EnzCap-minPro-ProwasteLim-prot_mass_fraction_max_eqn-membrane_SA_vivo-Weight_Prowaste-kappEstSlack-membrane_SA_eqn-v_trans_UB_fx-v_trans_LB_fx-v_trans_AA_sum_def-v_sum_def/;
rba.optfile = 1;
$offMulti

$if exist %phenotype_path% $include %phenotype_path%

alias(j_fva,j);

rxns_to_test(j)$blocked(j) = no;
* minimize dummy objective 1st, to save time optimizing
Solve rba using lp minimizing z_dummy;
rxns_to_test(j)$((v.l(j) gt %vmin% or v.l(j) lt -%vmin%) and (rxns_to_test(j))) = no;
*while(card(rxns_to_test) > 0,
loop(j_fva,
	put log; put 'Performing FVA for ' card(rxns_to_test) ' rxns'; putclose;
	if(card(rxns_to_test) le 0,
		break;
	);
	Solve rba using lp maximizing z_fva;
	if(rba.modelstat ne 1,
		put log; put 'SOLVESTAT=' rba.solvestat:0:0; putclose;
		if(rba.modelstat eq 13,
			rba.solvelink = 1;
			solve rba using lp minimizing z_dummy;
			if(rba.modelstat eq 1,
				rxns_to_test(j)$((v.l(j) gt %vmin% or v.l(j) lt -%vmin%) and (rxns_to_test(j))) = no;
				continue;
			);
		);
		break;
	);
	if((sum(j$((v.l(j) gt %vmin% or v.l(j) lt -%vmin%) and (rxns_to_test(j))),1) eq 0),
		break;
	);
	rxns_to_test(j)$((v.l(j) gt %vmin% or v.l(j) lt -%vmin%) and (rxns_to_test(j))) = no;
);
rxns_to_test(j)$((v.l(j) gt %vmin% or v.l(j) lt -%vmin%) and (rxns_to_test(j))) = no;
* test out solving w/o these rxns
v.fx(j)$rxns_to_test(j)=0;
solve rba using lp minimizing z_dummy;
file ffva /rxns_to_test.txt/;
put ffva;
put '/'/;
loop(j$rxns_to_test(j),
	put "'" j.tl:0 "'"/;
);
put '/';
if (rba.modelStat ne 1, abort.noError "blocking remaining fluxes prevents an optimal solution from being found";);

file fblock /blocked_rxns.txt/;
put fblock;
put '/'/;
loop(j$rxns_to_test(j),
	put "'" j.tl:0 "'"/;
);
put '/';

$include %model_root_path%GAMS/default_outputs.gms