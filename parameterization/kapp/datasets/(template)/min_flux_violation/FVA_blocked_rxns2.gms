* Run FVA (Flux Variability Analysis) on the RBA model, using the RBA model's structure to prioritize which reactions to test first and speeding up the process.
* This does involve blocking rxns whose enzymes can't be allocated to them.
* To use in conjunction with another script, set the objective function to its optimal value from that script and then include this file
$if not setGlobal gms $setGlobal gms %system.FN%

$setGlobal rerun 0
$if exist "rxns_to_test.txt" $setGlobal rerun 1

$INLINECOM /*  */
$include ../paths.txt
$offInline

* run RBA without bisection, using default settings and outputs
$include kapp_calc_defaults.gms

*prosynSlackUB.up(pro_measured) = inf;
*prosynSlackLB.up(pro_measured) = 1 + %bi_tol%;

$onMultiR
Sets
rxns_to_test(j)
$ifThen %rerun%==1
$include "rxns_to_test.txt"
$else $include %rxns_path%
$endIf
priority_stage /prosyn, ptm, enzsyn, metabolism, other/
rxns_priority(j,priority_stage)
rxns_to_prioritize(j) /#j/
i_pro(i)
$include %gms_path%pro.txt
i_ptmpro(i)
$include %gms_path%ptmpro.txt
blocked(j) /
$ifThen %ignore_blocked_j%==0
$if exist "blocked_rxns.txt" $include "blocked_rxns.txt"
$endIf
/
not_blocked(j) //
;
$offMulti

rxns_priority(j,'prosyn') = prosyn(j);
rxns_priority(j,'ptm') = ptm(j);
rxns_priority(j,'enzsyn') = rxns_enzsyn(j);
rxns_priority(j,'metabolism') = rxns_metab(j);
* add rxns that aren't found in the other categories to the 'other' category, in a way that automatically adapts to new priorities being added to the model
rxns_priority(j,'other') = not (rxns_priority(j,'prosyn') or rxns_priority(j,'ptm') or rxns_priority(j,'enzsyn') or rxns_priority(j,'metabolism'));

loop(enzload_rxn_coupling(j,j1),
    if (v.up(j1) eq 0,
        v.up(j) = 0;
    );
);

$include %model_root_path%find_blocked_rxns.gms

variables z_fva;

equation obj_fva; obj_fva.. z_fva =e= sum(j$rxns_to_prioritize(j),v(j));
$onMultiR
model rba /all-EnzCap-minPro-ProwasteLim-prot_mass_fraction_max_eqn-membrane_SA_vivo-Weight_Prowaste-kappEstSlack-membrane_SA_eqn-v_trans_UB_fx-v_trans_LB_fx-v_trans_AA_sum_def-v_sum_def/;
rba.optfile = 1;
$offMulti

$if exist %phenotype_path% $include %phenotype_path%

alias(j_fva,j);

rxns_to_test(j)$((v.lo(j) gt -%vmin% and v.up(j) lt %vmin%) or (v.lo(j) eq v.up(j)) or (sum(sm_j$dir(sm_j,j),1) ge 1 and sum(sm_j$dir(sm_j,j),v_sm_lo(sm_j)-v_sm_up(sm_j)) eq 0)) = no;
* minimize dummy objective 1st, to save time optimizing
Solve rba using lp minimizing z_dummy;
rxns_to_test(j)$((v.l(j) gt %vmin% or v.l(j) lt -%vmin%) and (rxns_to_test(j))) = no;
loop(priority_stage,
	rxns_to_prioritize(j) = rxns_priority(j,priority_stage) and rxns_to_test(j);
	if(card(rxns_to_prioritize) le 0,
		put log; put 'No rxns to test for priority stage ' priority_stage.tl:0; putclose;
		continue;
	);
	loop(j_fva,
		rxns_to_prioritize(j)$(not rxns_to_test(j)) = no;
		put log; put 'Performing FVA for ' priority_stage.tl:0 ' (' card(rxns_to_test) ' rxns remaining)'; putclose;
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
		if((sum(j$((v.l(j) gt %vmin% or v.l(j) lt -%vmin%) and rxns_to_prioritize(j)),1) eq 0),
* block all rxns dependent on this stage (assumes that only 1 rxn makes a given protein or post-translationally modified protein)
			blocked(j)$(rxns_to_prioritize(j)) = yes;
			if(sameas(priority_stage,'prosyn'),
				loop(i_pro,
					if(sum(j$(prosyn(j) and blocked(j) and S(i_pro,j) gt 0),1) gt 0,
						blocked(j)$(S(i_pro,j) lt 0) = yes;
					);
				);
			);
			if(sameas(priority_stage,'prosyn') or sameas(priority_stage,'ptm'),
				loop(i_ptmpro,
					if(sum(j$(ptm(j) and blocked(j) and S(i_ptmpro,j) gt 0),1) gt 0,
						blocked(j)$(S(i_ptmpro,j) lt 0) = yes;
					);
				);
				blocked(j)$(sum(j1$(blocked(j1) and pro_syn_waste(j1,j)),1) eq 1) = yes;
			);
			if(not sameas(priority_stage,'metabolism'),
				loop(i_enz,
					if(sum(j$(not blocked(j) and S(i_enz,j) gt 0),1) eq 0,
						blocked(j)$(S(i_enz,j) ne 0) = yes;
					);
				);
				loop(j$(rxns_with_enz(j) and not blocked(j)),
					if(sum(j1$(not blocked(j1) and enzload_rxn_coupling(j1,j)),1) eq 0,
						blocked(j) = yes;
					);
				);
			);
			rxns_to_test(j)$blocked(j) = no;
			break;
		);
		rxns_to_test(j)$((v.l(j) gt %vmin% or v.l(j) lt -%vmin%) and (rxns_to_test(j))) = no;
	);
);
blocked(j)$(sum(j1$(blocked(j1) and pro_syn_waste(j1,j)),1) eq 1) = yes;
rxns_to_test(j)$(((v.l(j) gt %vmin% or v.l(j) lt -%vmin%) or blocked(j)) and (rxns_to_test(j))) = no;
* test out solving w/o these rxns
parameter rxns_to_test_UB(j) "upper bounds for rxns to test";
rxns_to_test_UB(j) = v.up(j);
v.fx(j)$rxns_to_test(j)=0;
solve rba using lp minimizing z_dummy;
file ffva /rxns_to_test.txt/;
put ffva;
put '/'/;
loop(j$rxns_to_test(j),
	put "'" j.tl:0 "'"/;
);
put '/';
putclose;

file fblock /blocked_rxns.txt/;
put fblock;
$ifThen %rerun%==1
fblock.ap=1;
$endIf
loop(j$(blocked(j)),
	put "'" j.tl:0 "'"/;
);
if (rba.modelStat ne 1, abort.noError "blocking remaining fluxes prevents an optimal solution from being found";);

* do FVA on the remaining rxns_to_test, to see if they are truly blocked
v.up(j)$rxns_to_test(j) = rxns_to_test_UB(j);
loop(j_fva$rxns_to_test(j_fva),
	if(not_blocked(j_fva),
		continue;
	);
	rxns_to_prioritize(j) = no;
	rxns_to_prioritize(j_fva) = yes;
    Solve rba using lp maximizing z_fva;
    if(rba.modelstat ne 1,
        put log; put 'SOLVESTAT=' rba.solvestat:0:0; putclose;
        abort.noError "FVA on remaining rxns failed";
    );
	if((v.l(j_fva) le %vmin% and v.l(j_fva) ge -%vmin%),
		put "'" j_fva.tl:0 "'"/;
		v.fx(j_fva) = 0;
	);
	not_blocked(j)$(rxns_to_test(j) and v.l(j) gt %vmin%) = yes;
);
* more sequential blocks
loop(i_pro,
	if(sum(j$(prosyn(j) and blocked(j) and S(i_pro,j) gt 0),1) gt 0,
		loop(j$(S(i_pro,j) lt 0 and not blocked(j)),
			put "'" j.tl:0 "'"/;
			blocked(j) = yes;
		);
	);
);
loop(i_ptmpro,
	if(sum(j$(ptm(j) and blocked(j) and S(i_ptmpro,j) gt 0),1) gt 0,
		loop(j$(S(i_ptmpro,j) lt 0 and not blocked(j)),
			put "'" j.tl:0 "'"/;
			blocked(j) = yes;
		);
	);
);
loop(i_enz,
	if(sum(j$(not blocked(j) and S(i_enz,j) gt 0),1) eq 0,
		loop(j$(S(i_enz,j) ne 0 and not blocked(j)),
			put "'" j.tl:0 "'"/;
			blocked(j) = yes;
		);
	);
);
loop(j$(rxns_with_enz(j) and not blocked(j)),
	if(sum(j1$(not blocked(j1) and enzload_rxn_coupling(j1,j)),1) eq 0,
		put "'" j.tl:0 "'"/;
		blocked(j) = yes;
	);
);

$include %model_root_path%GAMS/default_outputs.gms

