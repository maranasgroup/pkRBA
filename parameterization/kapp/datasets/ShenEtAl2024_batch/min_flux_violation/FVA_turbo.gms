* Run FVA (Flux Variability Analysis) on the RBA model, using the RBA model's structure to prioritize which reactions to test first and speeding up the process.
* This does involve blocking rxns whose enzymes can't be allocated to them.
* To use in conjunction with another script, set the objective function to its optimal value from that script and then include this file
$if not setGlobal gms $setGlobal gms %system.FN%

* turn off unnecessary rxns to improve solving time
$ifThen %gms%==FVA_turbo
$INLINECOM /*  */
$include ../paths.txt
$offInline

* run RBA without bisection, using default settings and outputs
$include kapp_calc_defaults.gms

v.fx(j)$enzsyn(j) = 0;
v.fx(j)$enzload(j) = 0;
$endIf

Sets
rxns_to_test(j) "rxns to include in FVA" /#j
/
priority_stage /prosyn, ptm, enzsyn, metabolism, misc/
rxns_priority(j,priority_stage)
rxns_to_prioritize(j) /#j/
dummy(j) "dummy set for iterating over rxns"
i_pro(i)
$include %gms_path%pro.txt
i_ptmpro(i)
$include %gms_path%ptmpro.txt
;

rxns_priority(j,'prosyn') = prosyn(j);
rxns_priority(j,'ptm') = ptm(j);
rxns_priority(j,'enzsyn') = rxns_enzsyn(j);
rxns_priority(j,'metabolism') = rxns_metab(j);
*rxns_priority(j,'misc') = yes$(sum(priority_stage$(rxns_priority(j,priority_stage)),1) eq 0);
rxns_priority(j,'misc') = not (rxns_priority(j,'prosyn') or rxns_priority(j,'ptm') or rxns_priority(j,'enzsyn') or rxns_priority(j,'metabolism'));
put log; put 'misc set'/; putclose;
loop(j$rxns_priority(j,'misc'),
	put log; put j.tl:0/; putclose;
);

loop(enzload_rxn_coupling(j,j1),
	if (v.up(j1) eq 0,
		v.up(j) = 0;
	);
);

Parameters
high(j), low(j), high_stat(j), low_stat(j)
;

Variables
z_fva
;

Equations
Obj_FVA
;

Obj_FVA..           z_fva =e= sum(j$dummy(j),v(j));
$onMultiR
*model rba /all-ProsynSlacks-sm_LB_exp-sm_UB_exp-fluxSlackBounds-lumped_flux-FlexibleRiboCapacity-Obj1-Obj2-Obj3-EnzCap-minPro-ProwasteLim-prot_mass_fraction_max_eqn-membrane_SA_vivo-Weight_Prowaste-kappEstSlack-Obj5-membrane_SA_eqn-v_trans_UB_fx-v_trans_LB_fx-v_trans_AA_sum_def-v_sum_def/;
*model rba /all-ProsynSlacks-EnzCap-minPro-ProwasteLim-prot_mass_fraction_max_eqn-membrane_SA_vivo-Weight_Prowaste-kappEstSlack-Obj5-membrane_SA_eqn-v_trans_UB_fx-v_trans_LB_fx-v_trans_AA_sum_def-v_sum_def/;
model rba /all-EnzCap-minPro-ProwasteLim-prot_mass_fraction_max_eqn-membrane_SA_vivo-Weight_Prowaste-membrane_SA_eqn-v_trans_UB_fx-v_trans_LB_fx-v_trans_AA_sum_def-v_sum_def/;
rba.optfile = 1;
$offMulti

*Updates output while solving, so errors won't get rid of output
file f_fva /%gms%.txt/;
put f_fva;
* Clear old file
put "rxn":66, "min_flux":15, "max_flux":15, "high_stat":10, "low_stat":10/;
* Append old file
f_fva.ap=1;
f_fva.nr=2;

$if exist %phenotype_path% $include %phenotype_path%

v.up(j)$(sum(i_pro$(sum(j1$(prosyn(j1) and v.up(j1) lt %vmin% and S(i_pro,j1) gt 0),1) gt 0 and S(i_pro,j) lt 0),1) gt 0 and not (v.up(j) lt %vmin%)) = 0;
v.up(j)$(sum(i_ptmpro$(sum(j1$(ptm(j1) and v.up(j1) lt %vmin% and S(i_ptmpro,j1) gt 0),1) gt 0 and S(i_ptmpro,j) lt 0),1) gt 0 and not (v.up(j) lt %vmin%)) = 0;
v.up(j)$(sum(j1$(v.up(j1) lt %vmin% and pro_syn_waste(j1,j)),1) eq 1) = v.lo(j);
v.up(j)$(sum(i_enz$(sum(j1$(not (v.up(j1) lt %vmin%) and S(i_enz,j1) gt 0),1) eq 0 and S(i_enz,j) ne 0),1) gt 0 and not (v.up(j) lt %vmin%)) = 0;
v.up(j)$(rxns_with_enz(j) and not (v.up(j) lt %vmin%) and sum(j1$(not (v.up(j1) lt %vmin%) and enzload_rxn_coupling(j1,j)),1) eq 0) = 0;

rxns_to_test(j)$((v.lo(j) gt -%vmin% and v.up(j) lt %vmin%) or (v.lo(j) eq v.up(j))) = no;
alias(j2,j);

loop(priority_stage,
	rxns_to_prioritize(j) = rxns_priority(j,priority_stage) and rxns_to_test(j) and (v.up(j) gt %vmin%) and (v.up(j) ne v.lo(j));
	if(card(rxns_to_prioritize) le 0,
		put log; put 'No rxns to test for priority stage ' priority_stage.tl:0; putclose;
		continue;
	);
	loop(j2$(rxns_to_prioritize(j2)),
		if(v.lo(j2) ne v.up(j2),
			dummy(j) = no;
			dummy(j2) = yes;
			high(j2) = v.lo(j2);
			low(j2) = v.lo(j2);
			high_stat(j2) = -1;
			low_stat(j2) = -1;
			Solve rba using lp maximizing z_fva;
			high_stat(j2) = rba.modelstat;
			If((rba.modelstat eq 1),
				if(z_fva.l lt %vmin%,
					v.up(j2) = max(z_fva.l,0);
					v.fx(j)$(sum(i_pro$(sum(j1$(prosyn(j1) and v.up(j1) lt %vmin% and S(i_pro,j1) gt 0),1) gt 0 and S(i_pro,j) lt 0),1) gt 0 and not (v.up(j) lt %vmin%)) = 0;
					v.fx(j)$(sum(i_ptmpro$(sum(j1$(ptm(j1) and v.up(j1) lt %vmin% and S(i_ptmpro,j1) gt 0),1) gt 0 and S(i_ptmpro,j) lt 0),1) gt 0 and not (v.up(j) lt %vmin%)) = 0;
					v.up(j)$(sum(j1$(v.up(j1) lt %vmin% and pro_syn_waste(j1,j)),1) eq 1) = v.lo(j);
					v.fx(j)$(sum(i_enz$(sum(j1$(not (v.up(j1) lt %vmin%) and S(i_enz,j1) gt 0),1) eq 0 and S(i_enz,j) ne 0),1) gt 0 and not (v.up(j) lt %vmin%)) = 0;
					v.fx(j)$(rxns_with_enz(j) and not (v.up(j) lt %vmin%) and sum(j1$(not (v.up(j1) lt %vmin%) and enzload_rxn_coupling(j1,j)),1) eq 0) = 0;
				);
				if((z_fva.l gt v.lo(j2)),
					high(j2) = z_fva.l;
					solve rba using lp minimizing z_fva;
					low_stat(j2) = rba.modelstat;
					if(rba.modelstat eq 1,
						low(j2) = z_fva.l;
					);
				else
					low_stat(j2) = 1;
				);
			);
			put f_fva;
			put j2.tl:66, (low(j2) / %nscale%):15:8, (high(j2) / %nscale%):15:8, high_stat(j2):10:0, low_stat(j2):10:0/;
			putclose f_fva;
		);
	);
);

file f_fva2 /%gms%.bounds.txt/;
f_fva2.ap=1;
f_fva2.nr=2;
put f_fva2;
put "":66, "min_flux":15, "max_flux":15/;
loop(j2$(rxns_to_test(j2)),
	if(((v.lo(j2) ne 0) or (v.up(j2) ne 0)) and (high_stat(j2) eq 1) and (low_stat(j2) eq 1),
		put j2.tl:66, low(j2):15:8, high(j2):15:8/;
	);
);
putclose f_fva2;

* separate these to distinguish tested rxns from those that were blocked (i.e., not tested)
file f_fva3 /%gms%.bounds_all.txt/;
f_fva3.nr=2;
put f_fva3;
put "":66, "lo":15, "up":15/;
v.up(j)$(high_stat(j) eq 1) = high(j);
v.lo(j)$(low_stat(j) eq 1) = low(j);
loop(j2,
	put j2.tl:66, v.lo(j2):15:8, v.up(j2):15:8/;
);
putclose f_fva3;

$include %model_root_path%GAMS/default_outputs.gms
