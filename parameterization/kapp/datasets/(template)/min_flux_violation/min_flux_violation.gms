* Minimize violation of fluxes that are zero due to non-production of enzymes

* define current gms file for use in other files
$setGlobal gms %system.FN%

* load default sets, parameters, and equations for kapp calculation
$include kapp_calc_defaults.gms

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

Parameters
prosynSlackSum_opt
fluxSlack_opt
fluxSum_j_NP_opt
fluxSum_opt
kappEstSlackSum_opt
w_prowaste_opt
;

$include %model_root_path%find_blocked_rxns.gms

* minimize disagreement with proteomics data, while allowing some where needed (e.g., measurement errors)
* exclude v_trans equations since they're not needed and may slow down run times; add back if using them here
Model kapp_calc /all-EnzCap-minPro-ProwasteLim-prot_mass_fraction_max_eqn-membrane_SA_vivo-Weight_Prowaste-kappEstSlack-membrane_SA_eqn-v_trans_UB_fx-v_trans_LB_fx-v_trans_AA_sum_def-v_sum_def/;
kapp_calc.optfile = 1;
Solve kapp_calc using lp minimizing prosynSlackSum;
put log; put 'minimized prosynSlackSum'/; putclose;

* initially allow some minor slack (to account for measurement error) on just proteins w/o modeled functions, to minimize the potential impact on enzymes
* do positive 1st, to minimize impact on ribosome calculations
if (kapp_calc.modelStat ne 1, 
    put log; put 'Model infeasible - relaxing proteomics constraints'/; putclose;
	prosynSlackUB.up(pro_measured)$(pro_wmf(pro_measured) and not unlimited_slack(pro_measured)) = 0.05;
    Solve kapp_calc using lp minimizing prosynSlackSum;
);
* Try max possible positive slack for proteins w/o modeled functions, to avoid excluding proteins from consideration as ribo subunits
if (kapp_calc.modelStat ne 1, 
	prosynSlackUB.up(pro_measured)$pro_wmf(pro_measured) = inf;
	Solve kapp_calc using lp minimizing prosynSlackSum;
);
* Try adding small slack to all proteins, to account for measurement error
if (kapp_calc.modelStat ne 1, 
	prosynSlackUB.up(pro_measured)$((not pro_wmf(pro_measured)) and (not unlimited_slack(pro_measured))) = 0.05;
	Solve kapp_calc using lp minimizing prosynSlackSum;
);
* Allow max positive slack to all proteins
if (kapp_calc.modelStat ne 1, 
	prosynSlackUB.up(pro_measured) = inf;
	Solve kapp_calc using lp minimizing prosynSlackSum;
);
* Allow negative slack
if (kapp_calc.modelStat ne 1, 
* set bounds to min. height needed to cover full range w/o being excessive (to reduce solving time)
	prosynSlackLB.up(pro_measured) = 1 + %bi_tol%;
	Solve kapp_calc using lp minimizing prosynSlackSum;
);
if (kapp_calc.modelStat ne 1, put ff; put kapp_calc.modelStat/; putclose; abort.noError "no optimal solution found even after relaxing constraints";);

prosynSlackSum_opt = prosynSlackSum.l;
prosynSlackSum.up = max(prosynSlackSum_opt * (1+%epsilon%), prosynSlackSum_opt + %epsilon%);
* force each variable in sum to not exceed current sum (+ some leeway for numerical precision) or its prior bounds (to prevent this from adding flexibility where unexpected), to improve solving times and reduce degeneracy
prosynSlackUB.up(pro_measured) = min(prosynSlackSum_opt + (2 * max(%vmin%,%epsilon%)), prosynSlackUB.up(pro_measured));
prosynSlackLB.up(pro_measured) = min(prosynSlackSum_opt + (2 * max(%vmin%,%epsilon%)), prosynSlackLB.up(pro_measured));

Solve kapp_calc using lp minimizing fluxSlack;
put log; put 'minimized fluxSlack'/; putclose;
if (kapp_calc.modelStat ne 1, abort.noError "no optimal solution found";);
fluxSlack_opt = fluxSlack.l;
fluxSlack.up = fluxSlack_opt * (1+%epsilon%);
* force each variable in sum to not exceed current sum (+ some leeway for numerical precision) or its prior bounds (to prevent this from adding flexibility where unexpected), to improve solving times and reduce degeneracy
s_v_exp_ub.up(sm_j_measured) = min(fluxSlack.l + (2 * max(%vmin%,%epsilon%)), s_v_exp_ub.up(sm_j_measured));
s_v_exp_lb.up(sm_j_measured) = min(fluxSlack.l + (2 * max(%vmin%,%epsilon%)), s_v_exp_lb.up(sm_j_measured));

Solve kapp_calc using lp minimizing fluxSum_j_NP;
put log; put 'minimized fluxSum_j_NP'/; putclose;
if (kapp_calc.modelStat ne 1, abort.noError "no optimal solution found";);
* force rxns that were turned off to stay off
*v.fx(j)$(rxns_with_no_prodata(j) and (v.l(j) eq 0)) = 0;
fluxSum_j_NP_opt = fluxSum_j_NP.l;
fluxSum_j_NP.up = fluxSum_j_NP_opt * (1+%epsilon%);

* Put after minimizing unmeasured fluxes since not all PTMs may be crucial for growth, but we want to capture as many as possible
* Solve kapp_calc using lp maximizing SumPTM;
* put log; put 'maximized SumPTM'/; putclose;
* if (kapp_calc.modelStat ne 1, abort.noError "no optimal solution found";);
* SumPTM.lo = SumPTM.l - %epsilon%;

Solve kapp_calc using lp minimizing fluxSum;
put log; put 'minimized fluxSum'/; putclose;
if (kapp_calc.modelStat ne 1, abort.noError "no optimal solution found";);
* force total flux to previous value
fluxSum_opt = fluxSum.l;
fluxSum.up = fluxSum_opt*(1+%tol%);

* Solve again, encouraging more equal use of all enzymes
* NOTE: disabled this step since it can lead to arbitrary reductions in ENZLOAD fluxes, even when other ones aren't being used. This can increase kapps by reducing the denominator; how to fix this is unclear.
*Solve kapp_calc using lp minimizing slackSum;
*put log; put 'minimized uneven enzload distribution'/; putclose;

ff.nr = 2; put ff;
put kapp_calc.modelStat/;
putclose;

file ff2 /%system.FN%.objectives.txt/;
ff2.nr = 2; put ff2; ff2.pc=6;
put 'prosynSlackSum',prosynSlackSum.l:0:11/;
put 'fluxSlackSum',(fluxSlack.l/%nscale%):0:11/;
put 'fluxSum_j_NP',(fluxSum_j_NP.l/%nscale%):0:11/;
put 'fluxSum',(fluxSum.l/%nscale%):0:11/;
putclose;

file ff3 /%system.FN%.flux_gamsscaled.txt/;
ff3.nr = 2; put ff3;
loop(j,
    if ( (v.l(j) gt %vmin%),
        put j.tl:0, system.tab, 'v', system.tab, v.l(j):0:15/;
    );
);
putclose;

file ff3a /%system.FN%.flux_unscaled.txt/;
ff3a.nr = 2; put ff3a;
loop(j,
    if ( (v.l(j) gt %vmin%),
        put j.tl:0, system.tab, 'v', system.tab, (v.l(j)/%nscale%):0:15/;
    );
);
putclose;

file ff4 /%system.FN%.flux_essential_with_no_prodata_gamsscaled.txt/;
file ff4c /%system.FN%.rxns_essential_with_no_prodata_gamsscaled.txt/;
ff4.nr = 2; put ff4;
ff4c.nr = 2; put ff4c;
put ff4c '/'/;
loop(j$rxns_with_no_prodata(j),
    if ( (v.l(j) gt %vmin%),
        put ff4 j.tl:0, system.tab, 'v', system.tab, v.l(j):0:15/ ff4c "'" j.tl:0 "'" system.tab v.l(j):0:15/;
    );
);
put ff4c '/';
putclose ff4 ff4c;

file ff4a /%system.FN%.flux_essential_with_no_prodata_unscaled.txt/;
ff4a.nr = 2; put ff4a;
loop(j$rxns_with_no_prodata(j),
    if ( (v.l(j) gt %vmin%),
        put j.tl:0, system.tab, 'v', system.tab, (v.l(j)/%nscale%):0:15/;
    );
);
putclose;

file ff4b /%system.FN%.rxns_nonessential_with_no_prodata.txt/;
ff4b.nr = 2; put ff4b;
put '/'/;
loop(j$(rxns_with_no_prodata(j)),
    if ( (v.l(j) le %vmin%),
        put j.tl:0/;
    );
);
put '/'/;
putclose;

file ff6 /%system.FN%.s_v_exp.txt/;
ff6.nr = 2; ff6.pc=6; put ff6;
loop(sm_j_measured,
    if ( (s_v_exp_lb.l(sm_j_measured) gt %vmin%) or (s_v_exp_ub.l(sm_j_measured) gt %vmin%),
        put sm_j_measured.tl:0, s_v_exp_lb.l(sm_j_measured):0:15, s_v_exp_ub.l(sm_j_measured):0:15/;
    );
);
loop(sm_j_lumped,
    if ( (s_v_exp_lb_lumped.l(sm_j_lumped) gt %vmin%) or (s_v_exp_ub_lumped.l(sm_j_lumped) gt %vmin%),
        put sm_j_lumped.tl:0, s_v_exp_lb_lumped.l(sm_j_lumped):0:15, s_v_exp_ub_lumped.l(sm_j_lumped):0:15/;
    );
);
putclose;

file ff7 /%system.FN%.prosynSlack.txt/;
ff7.nr = 2; put ff7;
put 'index', system.tab, 'prosynSlack(% higher/lower than measured value)'/; 
loop(pro_measured,
    if ( prosynSlackUB.l(pro_measured) > %vmin% or prosynSlackLB.l(pro_measured) > %vmin%,
        put pro_measured.tl:0, system.tab, (100*prosynSlackUB.l(pro_measured)-100*prosynSlackLB.l(pro_measured)):0:15/;
    );
);
putclose;

file enz_flux /%system.FN%.enz_flux_calculation.txt/;
enz_flux.nr = 2; put enz_flux;
loop(j$(rxns_enzsyn(j) or rxns_enzload(j)),
    if ( (v.l(j) gt %vmin%),
        put j.tl:0, system.tab, (v.l(j)/%nscale%):0:18/;
    );
);
putclose;

* output protein levels for enforcing kapp calculation levels when needed
file ff9 /%system.FN%.prosyn_unscaled.txt/;
ff9.nr=2; ff9.nz=1e-30; put ff9;
put '/'/;
loop(j$prosyn(j),
    put "'" j.tl:0 "' " (v.l(j)/%nscale%):0:15/;
);
put '/'/;
putclose;

* output protein levels for enforcing kapp calculation levels when needed
file ff10 /%system.FN%.prosyn_gamsscaled.txt/;
ff10.nr=2; ff10.nz=1e-30; put ff10;
put '/'/;
loop(j$prosyn(j),
    put "'" j.tl:0 "' " v.l(j):0:15/;
);
put '/'/;
putclose;

* output protein levels for enforcing kapp calculation levels when needed
file ff11 /%system.FN%.prosyn_frac_unscaled.txt/;
ff11.nr=2; ff11.nz=1e-30; put ff11;
put '/'/;
loop(j$prosyn(j),
    put "'" j.tl:0 "' " (v.l(j)/v.l('BIOSYN-PROTTOBIO')):0:15/;
);
put '/'/;
putclose;

* output protein levels for enforcing kapp calculation levels when needed
file ff12 /%system.FN%.prosyn_nonzero_gamsscaled.txt/;
ff12.nr=2; ff12.nz=1e-30; put ff12;
put '/'/;
loop(j$prosyn(j),
	if ( (v.l(j) ge %vmin%),	
		put "'" j.tl:0 "' " v.l(j):0:15/;
	);
);
put '/'/;
putclose;
