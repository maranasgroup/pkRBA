* identical to enz_alloc.gms but meant to test isozymes
* finalize enzyme allocation

* define current gms file for use in other files
$if not setGlobal gms $setGlobal gms %system.FN%

* load default sets, parameters, and equations for kapp calculation
$include kapp_calc_defaults.gms

Sets
rxns_enzload_used(j)
$include "%rxns_enzload_used_path%"
rxns_NP_nonessential(j)
$include "./min_flux_violation.rxns_nonessential_with_no_prodata.txt"
objective_fluxes(j) / 'BIOSYN-PROTTOBIO' /
rxns_enzload_nonused(j) "enzload reactions for enzymatic rxns that aren't being used"
;
rxns_enzload_nonused(j) = rxns_enzload(j) and not rxns_enzload_used(j);

Parameters
kapp(j)
$include "./kapp_init.txt"
kappMaxSlackSum_opt
prosynSlackSum_opt
fluxSlack_opt
fluxSum_j_NP_opt
fluxSum_opt
kappEstSlackSum_opt
w_prowaste_opt
;
Variables
fluxSum_enzload_nonused
kappMaxSlack(rxns_with_enz)
kappMaxSlackSum
enzload_is_used(enzload)
enzload_is_used_sum
enzload_sum
;
kappMaxSlack.lo(rxns_with_enz) = 0;
kappMaxSlack.up(rxns_with_enz) = inf;
kappMaxSlackSum.lo = 0;
enzload_is_used.up(enzload) = 1;
enzload_is_used_sum.up = card(enzload);
enzload_sum.up = sum(enzload, v.up(enzload));
$ifThen setGlobal kappMaxSlackSum_Max
kappMaxSlackSum.up = %kappMaxSlackSum_Max%;
kappMaxSlack.up(rxns_with_enz) = min(%kappMaxSlackSum_Max%,kappMaxSlack.up(rxns_with_enz));
$endIf

loop(enzload_rxn_coupling(j,j1),
    if (v.up(j1) eq 0,
        v.up(j) = 0;
    );
);
v.up(j)$(rxns_with_enz(j) and sum(j1$(enzload_rxn_coupling(j1,j) and v.up(j1) gt 0),1) eq 0) = 0;
$include %model_root_path%find_blocked_rxns.gms
loop(enzload_rxn_coupling(j,j1),
    if (v.up(j1) eq 0,
        v.up(j) = 0;
    );
);
v.up(j)$(rxns_with_enz(j) and sum(j1$(enzload_rxn_coupling(j1,j) and v.up(j1) gt 0),1) eq 0) = 0;

* slacks for apportioning enzyme load according to flux, to account for unmeasured enzymes that may not have been produced (causes errors, so consider avoiding)
kappEstSlackPos.fx(rxns_with_enz) = 0;
kappEstSlackNeg.fx(rxns_with_enz) = 0;
*kappEstSlackPos.up(rxns_with_enz) = (v.up(rxns_with_enz) * mu) - %vmin%;
*kappEstSlackPos.up(rxns_with_enz)$(rxns_with_no_prodata(rxns_with_enz)) = sum(j1$enzload_rxn_coupling(j1,rxns_with_enz),%vmax% * %nscale% * kapp(j1)) + %epsilon%;
kappEstSlackPos.up(rxns_with_enz)$(v.up(rxns_with_enz) gt %vmin%) = sum(j1$enzload_rxn_coupling(j1,rxns_with_enz),v.up(j1) * kapp(j1)) - (%vmin% * mu) + %epsilon%;
kappEstSlackNeg.up(rxns_with_enz)$(v.up(rxns_with_enz) gt %vmin%) = max(v.up(rxns_with_enz) * mu - (smin(j1$(enzload_rxn_coupling(j1,rxns_with_enz) and v.up(j1) > %vmin%),kapp(j1))*%vmin%) + %epsilon%,0);
*kappEstSlackPos.lo(rxns_with_enz) = 0; kappEstSlackPos.up(rxns_with_enz) = inf;
*kappEstSlackNeg.lo(rxns_with_enz) = 0; kappEstSlackNeg.up(rxns_with_enz) = inf;
kappEstSlackPos.up(rxns_with_enz)$(v.up(rxns_with_enz) le %vmin%) = 0;
kappEstSlackNeg.up(rxns_with_enz)$(v.up(rxns_with_enz) le %vmin%) = 0;

** Disable enzyme load network for reactions that can't be used
*v.fx(j)$rxns_enzload(j) = 0;
* force production of all enzymes involved in used rxns
*v.lo(j)$rxns_enzload_used(j) = %vmin% * (1+1e-8); 
*v.up(j)$rxns_enzload_used(j) = %vmax% * %nscale%;

Equations
flux_enzload_nonused, kappEst(rxns_with_enz), kappMaxEst(rxns_with_enz), kappMaxSlackSum_def, enzload_sum_def
*, enzload_used_test(enzload), EIUsum_def
;

flux_enzload_nonused..         fluxSum_enzload_nonused =e= sum(j$rxns_enzload_nonused(j), v(j));
kappEst(rxns_with_enz)$(v.up(rxns_with_enz) gt %vmin%).. (v(rxns_with_enz) * mu) + kappEstSlackPos(rxns_with_enz) - kappEstSlackNeg(rxns_with_enz) =e= sum(j1$enzload_rxn_coupling(j1,rxns_with_enz),v(j1) * kapp(j1));
* max predicted kapp, for use in approximating enzyme levels: mu*vmax*nscale/vmin
kappMaxEst(rxns_with_enz)$(v.up(rxns_with_enz) gt %vmin% and kappMaxSlackSum.up gt %vmin%).. (v(rxns_with_enz) * mu) - kappMaxSlack(rxns_with_enz) =l= sum(j1$enzload_rxn_coupling(j1,rxns_with_enz), v(j1) * (mu * %vmax% * %nscale% / %vmin%));
kappMaxSlackSum_def..          kappMaxSlackSum =e= sum(rxns_with_enz, kappMaxSlack(rxns_with_enz));
enzload_sum_def..				enzload_sum =e= sum(j$enzload(j), v(j));
*enzload_used_test(enzload)..	enzload_is_used(enzload) =l= (v(enzload)-%vmin%)/%vmin%;
*EIUsum_def..	enzload_is_used_sum =e= sum(enzload, enzload_is_used(enzload));
*$include "./enz_alloc_constraints.txt"
* minimize disagreement with proteomics data, while allowing some where needed (e.g., measurement errors)
* exclude v_trans equations since they're not needed and may slow down run times; add back if using them here
Model m1 /all-EnzCap-minPro-ProwasteLim-prot_mass_fraction_max_eqn-membrane_SA_vivo-membrane_SA_eqn-v_trans_UB_fx-v_trans_LB_fx-v_trans_AA_sum_def/;
m1.optfile = 1;
m1.solveLink = 5;
if(kappMaxSlackSum.up gt %vmin%,
	Solve m1 using lp minimizing kappMaxSlackSum;

* initially allow some minor slack (to account for measurement error) on just proteins w/o modeled functions, to minimize the potential impact on enzymes
* do positive 1st, to minimize impact on ribosome calculations
	if (m1.modelStat ne 1,
		put log; put 'Model infeasible - relaxing proteomics constraints'/; putclose;
		prosynSlackUB.up(pro_measured)$(pro_wmf(pro_measured) and not unlimited_slack(pro_measured)) = 0.05;
		Solve m1 using lp minimizing kappMaxSlackSum;
	);
* Try max possible positive slack for proteins w/o modeled functions, to avoid excluding proteins from consideration as ribo subunits
	if (m1.modelStat ne 1,
		prosynSlackUB.up(pro_measured)$pro_wmf(pro_measured) = inf;
		Solve m1 using lp minimizing kappMaxSlackSum;
	);
* Try adding small slack to all proteins, to account for measurement error
	if (m1.modelStat ne 1,
		prosynSlackUB.up(pro_measured)$((not pro_wmf(pro_measured)) and (not unlimited_slack(pro_measured))) = 0.05;
		Solve m1 using lp minimizing kappMaxSlackSum;
	);
* Allow max positive slack to all proteins
	if (m1.modelStat ne 1,
		prosynSlackUB.up(pro_measured) = inf;
		Solve m1 using lp minimizing kappMaxSlackSum;
	);
* Allow negative slack
	if (m1.modelStat ne 1,
* set bounds to min. height needed to cover full range w/o being excessive (to reduce solving time)
		prosynSlackLB.up(pro_measured) = 1 + (10*%bi_tol%);
		Solve m1 using lp minimizing kappMaxSlackSum;
	);
	if (m1.modelStat ne 1, put ff; put m1.modelStat/; putclose; abort.noError "no optimal solution found even after relaxing constraints";);
	put log; put 'minimized kappMaxSlackSum'/; putclose;
	kappMaxSlackSum_opt = kappMaxSlackSum.l;
	kappMaxSlackSum.up = kappMaxSlackSum_opt * (1+%epsilon%);
*kappMaxSlack.up(rxns_with_enz) = kappMaxSlackSum.l * (1+%epsilon%);
);
file ff_kms /%system.FN%.kappMaxSlack.txt/;
ff_kms.nr=2; ff_kms.pc=6; put ff_kms;
put 'kappMaxSlackSum=' kappMaxSlackSum.l:0/ '.m' kappMaxSlackSum.m:0/;
put 'j','v','v_enzload','kapp_init','kapp est.','slack'/;
loop(enzload_rxn_coupling(j1,rxns_with_enz),
	if (kappMaxSlack.l(rxns_with_enz) gt 0,
		put rxns_with_enz.tl:0, v.l(rxns_with_enz):0:15, v.l(j1):0:15, kapp(rxns_with_enz):0, (v.l(rxns_with_enz)*mu/max(v.l(j1), 1e-20)):0:15, kappMaxSlack.l(rxns_with_enz)/;
	);
);
putclose;

if(smax(pro_measured, prosynSlackUB.up(pro_measured) + prosynSlackLB.up(pro_measured)) gt %vmin%,
	Solve m1 using lp minimizing prosynSlackSum;
	put log; put 'minimized prosynSlackSum'/; putclose;
	if (m1.modelStat ne 1, abort.noError "no optimal solution found";);
* add vmin * length of the variables covered by the objective function (+1 to account for objective variable itself), to avoid precision errors
	prosynSlackSum_opt = prosynSlackSum.l;
*prosynSlackSum.up = prosynSlackSum_opt + (%vmin%*(1+card(pro_measured)));
	prosynSlackSum.up = prosynSlackSum_opt + max(%epsilon%, prosynSlackSum_opt*%epsilon%);
* force each variable in sum to not exceed current sum (+ some leeway for numerical precision) or its prior bounds (to prevent this from adding flexibility where unexpected), to improve solving times and reduce degeneracy
	prosynSlackUB.up(pro_measured) = min(prosynSlackSum_opt + (2 * max(%vmin%,%epsilon%)), prosynSlackUB.up(pro_measured));
	prosynSlackLB.up(pro_measured) = min(prosynSlackSum_opt + (2 * max(%vmin%,%epsilon%)), prosynSlackLB.up(pro_measured));
*prosynSlackUB.up(pro_measured) = min(prosynSlackSum.l + (2 * max(%vmin%,%epsilon%)), prosynSlackUB.up(pro_measured));
*prosynSlackLB.up(pro_measured) = min(prosynSlackSum.l + (2 * max(%vmin%,%epsilon%)), prosynSlackLB.up(pro_measured));
);

equation m1opt; m1opt.. prosynSlackSum =l= prosynSlackSum_opt + max(%epsilon%, prosynSlackSum_opt*%epsilon%);

model m2 /m1, m1opt/;
m2.optfile=1;
m2.solvelink=5;

Solve m2 using lp minimizing fluxSlack;
put log; put 'minimized fluxSlack'/; putclose;
if (m2.modelStat ne 1, abort.noError "no optimal solution found";);
fluxSlack_opt = fluxSlack.l;
*fluxSlack.up = fluxSlack_opt * (1+%epsilon%) * card(sm_j_measured);
fluxSlack.up = fluxSlack_opt + (%epsilon% * card(sm_j_measured));
*fluxSlack.up = fluxSlack.l * (1+%epsilon%);
* force each variable in sum to not exceed current sum (+ some leeway for numerical precision) or its prior bounds (to prevent this from adding flexibility where unexpected), to improve solving times and reduce degeneracy
s_v_exp_ub.up(sm_j_measured) = min(fluxSlack.up, s_v_exp_ub.up(sm_j_measured));
s_v_exp_lb.up(sm_j_measured) = min(fluxSlack.up, s_v_exp_lb.up(sm_j_measured));

Solve m2 using lp minimizing w_prowaste;
put log; put 'minimized prowaste mass'/; putclose;
if (m2.modelStat ne 1, abort.noError "no optimal solution found";);
w_prowaste_opt = w_prowaste.l;
* force every pro/ptm waste to be within the optimal limit
v.up(j)$((prowaste(j) or ptmwaste(j)) and S('BIO-protwasted',j) gt 0) = min((w_prowaste_opt / S('BIO-protwasted',j)) + max(%epsilon%,%epsilon%*(w_prowaste_opt / S('BIO-protwasted',j))),v.up(j));
equation m2opt; m2opt.. w_prowaste =l= w_prowaste_opt + max(%epsilon%,w_prowaste_opt*%epsilon%);
