* Run FVA (Flux Variability Analysis) on the RBA model
* To use in conjunction with another script, set the objective function to its optimal value from that script and then include this file
$if not setGlobal gms $setGlobal gms %system.FN%

$setGlobal ignore_blocked 1
$setGlobal rerun 0

$INLINECOM /*  */
$include ../paths.txt
$offInline

* run RBA without bisection, using default settings and outputs
$include kapp_calc_defaults.gms

prosynSlackUB.up(pro_measured) = inf;
prosynSlackLB.up(pro_measured) = 1 + %bi_tol%;

* turn off unnecessary rxns to improve solving time
v.fx(j)$enzsyn(j) = 0;
v.fx(j)$enzload(j) = 0;

** find all rxns j with identical bounds, and turn all but 1 of them off to reduce solving time
loop(sm_j,
	if(sum(j$dir(sm_j,j), 1) gt 1,
		loop(j$((dir(sm_j,j) ne 0) and (v.up(j) ne 0)),
			if(sum(j1$(dir(sm_j,j1) eq dir(sm_j,j) and (v.up(j) eq v.up(j1)) and (v.lo(j) eq v.lo(j1))), 1) gt 1,
				if(rxns_with_no_prodata(j) or (sum(j1$((not rxns_with_no_prodata(j1)) and dir(sm_j,j1) eq dir(sm_j,j) and (v.up(j) eq v.up(j1)) and (v.lo(j) eq v.lo(j1))), 1) gt 1),
					v.fx(j) = 0;
				);
			);
		);
	);
);

$onMultiR
Sets
rxns_to_test(sm_j)
$ifThen %rerun%==1
$include "rxns_to_test_SM.txt"
$else $include %sm_rxns_path%
$endIf
blocked_SM(sm_j) //;
scalar nonzero_ub
$offMulti
opt_direction /1/;

loop(sm_j,
	if((v_sm_up(sm_j) le %vmin%) and (v_sm_lo(sm_j) ge -%vmin%),
		loop(j$(dir(sm_j,j) ne 0),
			v.fx(j) = 0;
		);
	else
		if(cap_gross_flux(sm_j) ne 0,
			if(v_sm_lo(sm_j) ge -%vmin%,
				v.fx(j)$(dir(sm_j,j) = -1) = 0;
			elseif(v_sm_up(sm_j) le %vmin%),
				v.fx(j)$(dir(sm_j,j) = 1) = 0;
			);
		);
		nonzero_ub = 0;
		loop(j$(dir(sm_j,j) ne 0),
			if(v.up(j) gt %vmin%,
				nonzero_ub = 1;
				break;
			);
		);
		if(nonzero_ub = 0,
			v_sm_lo(sm_j) = 0;
			v_sm_up(sm_j) = 0;
			put log; put 'blocked ' sm_j.tl:0 ' due to lack of v(j) entries with bounds over ' %vmin%; putclose;
		);
	);
);

$include %model_root_path%find_blocked_rxns.gms
blocked_SM(sm_j)$(v_sm_lo(sm_j) gt -%vmin% and v_sm_up(sm_j) lt %vmin%)=yes;

variables z_fva, z_fva2;

equation obj_fva; obj_fva.. z_fva =e= sum(sm_j$rxns_to_test(sm_j),sum(j$dir(sm_j,j),dir(sm_j,j)*v(j)));
equation obj_fva2; obj_fva2.. z_fva2 =e= sum(sm_j$(rxns_to_test(sm_j) and v_sm_up(sm_j) gt %vmin%),sum(j$dir(sm_j,j),dir(sm_j,j)*v(j))) - sum(sm_j$(rxns_to_test(sm_j) and v_sm_up(sm_j) le %vmin%),sum(j$dir(sm_j,j),dir(sm_j,j)*v(j)));
*equation obj_fva; obj_fva.. z_fva =e= sum(sm_j$rxns_to_test(sm_j),sum(j$dir(sm_j,j),v(j)));
$onMultiR
*model rba /all-ProsynSlacks-sm_LB_exp-sm_UB_exp-fluxSlackBounds-lumped_flux-FlexibleRiboCapacity-Obj1-Obj2-Obj3-EnzCap-minPro-ProwasteLim-prot_mass_fraction_max_eqn-membrane_SA_vivo-Weight_Prowaste-kappEstSlack-Obj5-membrane_SA_eqn-v_trans_UB_fx-v_trans_LB_fx-v_trans_AA_sum_def-v_sum_def/;
*model rba /all-ProsynSlacks-EnzCap-minPro-ProwasteLim-prot_mass_fraction_max_eqn-membrane_SA_vivo-Weight_Prowaste-kappEstSlack-Obj5-membrane_SA_eqn-v_trans_UB_fx-v_trans_LB_fx-v_trans_AA_sum_def-v_sum_def/;
model rba /all-EnzCap-minPro-ProwasteLim-prot_mass_fraction_max_eqn-membrane_SA_vivo-Weight_Prowaste-kappEstSlack-membrane_SA_eqn-v_trans_UB_fx-v_trans_LB_fx-v_trans_AA_sum_def-v_sum_def/;
rba.optfile = 1;
$offMulti

$if exist %phenotype_path% $include %phenotype_path%

alias(sm_fva,sm_j);

rxns_to_test(sm_j)$(blocked_SM(sm_j) or (v_sm_lo(sm_j) eq v_sm_up(sm_j))) = no;
*while(card(rxns_to_test) > 0,
loop(sm_fva,
	put log; put 'Performing FVA for ' card(rxns_to_test) ' rxns'; putclose;
	if(card(rxns_to_test) le 0,
		break;
	);
	if(opt_direction eq 1,
		Solve rba using lp maximizing z_fva;
		if(z_fva.l le %vmin% and z_fva.l ge -%vmin%,
			opt_direction=-1;
		);
	else
		solve rba using lp minimizing z_fva;
		if(z_fva.l le %vmin% and z_fva.l ge -%vmin% and (sum(sm_j$((sum(j$dir(sm_j,j),dir(sm_j,j)*v.l(j)) gt %vmin% or sum(j$dir(sm_j,j),dir(sm_j,j)*v.l(j)) lt -%vmin%) and (rxns_to_test(sm_j))),1) eq 0),
			solve rba using lp maximizing z_fva2;
			if(z_fva2.l le %vmin% and z_fva2.l ge -%vmin% and (sum(sm_j$((sum(j$dir(sm_j,j),dir(sm_j,j)*v.l(j)) gt %vmin% or sum(j$dir(sm_j,j),dir(sm_j,j)*v.l(j)) lt -%vmin%) and (rxns_to_test(sm_j))),1) eq 0),
				break;
			);
		);
	);
	if(rba.modelstat ne 1,
		put log; put 'SOLVESTAT=' rba.solvestat:0:0; putclose;
		if(rba.modelstat eq 13,
			rba.solvelink = 1;
			solve rba using lp minimizing z_dummy;
			if(rba.modelstat eq 1,
				rxns_to_test(sm_j)$((sum(j$dir(sm_j,j),dir(sm_j,j)*v.l(j)) gt %vmin% or sum(j$dir(sm_j,j),dir(sm_j,j)*v.l(j)) lt -%vmin%) and (rxns_to_test(sm_j))) = no;
				continue;
			);
		);
		break;
	);
	rxns_to_test(sm_j)$((sum(j$dir(sm_j,j),dir(sm_j,j)*v.l(j)) gt %vmin% or sum(j$dir(sm_j,j),dir(sm_j,j)*v.l(j)) lt -%vmin%) and (rxns_to_test(sm_j))) = no;
);
if(rba.modelstat eq 1,
	rxns_to_test(sm_j)$((sum(j$dir(sm_j,j),dir(sm_j,j)*v.l(j)) gt %vmin% or sum(j$dir(sm_j,j),dir(sm_j,j)*v.l(j)) lt -%vmin%) and (rxns_to_test(sm_j))) = no;
);
* test out solving w/o these rxns
v_sm_lo(sm_j)$rxns_to_test(sm_j)=0;
v_sm_up(sm_j)$rxns_to_test(sm_j)=0;
solve rba using lp minimizing z_dummy;
file ffva /rxns_to_test_SM.txt/;
put ffva;
put '/'/;
loop(sm_j$rxns_to_test(sm_j),
	put "'" sm_j.tl:0 "'"/;
);
put '/';
if (rba.modelStat ne 1, abort.noError "blocking remaining fluxes prevents an optimal solution from being found";);

file fblock /blocked_rxns_SM.txt/;
put fblock;
put '/'/;
loop(sm_j$rxns_to_test(sm_j),
	put "'" sm_j.tl:0 "'"/;
);
put '/';

$include %model_root_path%GAMS/default_outputs.gms

