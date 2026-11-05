* identical to enz_alloc.gms but meant to test isozymes
* finalize enzyme allocation

* define current gms file for use in other files
$setGlobal gms %system.FN%


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
$ifThen exist "enz_alloc.bounds.txt"
table bounds(j, *)
$if exist "enz_alloc.bounds.txt" $include "enz_alloc.bounds.txt"
;
v.lo(j)$(enzload(j)) = max((bounds(j, 'min_flux_max_enz') * %nscale%) * (1-max(1e-1,%vmin%)),v.lo(j));
v.up(j)$((enzload(j))) = max(min((bounds(j, 'max_flux_max_enz') * %nscale%) * (1+max(1e-1,%vmin%)),v.up(j)),v.lo(j));
$endIf

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
        v.fx(j) = 0;
    );
);
v.fx(j)$(rxns_with_enz(j) and sum(j1$(enzload_rxn_coupling(j1,j) and v.up(j1) gt 0),1) eq 0) = 0;
$include %model_root_path%find_blocked_rxns.gms
loop(enzload_rxn_coupling(j,j1),
    if (v.up(j1) eq 0,
        v.fx(j) = 0;
    );
);
v.fx(j)$(rxns_with_enz(j) and sum(j1$(enzload_rxn_coupling(j1,j) and v.up(j1) gt 0),1) eq 0) = 0;

* slacks for apportioning enzyme load according to flux, to account for unmeasured enzymes that may not have been produced (causes errors, so consider avoiding)
kappEstSlackPos.fx(rxns_with_enz) = 0;
kappEstSlackNeg.fx(rxns_with_enz) = 0;
*kappEstSlackPos.up(rxns_with_enz) = (v.up(rxns_with_enz) * mu) - %vmin%;
*kappEstSlackPos.up(rxns_with_enz)$(rxns_with_no_prodata(rxns_with_enz)) = sum(j1$enzload_rxn_coupling(j1,rxns_with_enz),%vmax% * %nscale% * kapp(j1)) + %epsilon%;
kappEstSlackPos.up(rxns_with_enz)$(v.up(rxns_with_enz) gt %vmin%) = max(sum(j1$enzload_rxn_coupling(j1,rxns_with_enz),v.up(j1) * kapp(j1)) - (%vmin% * mu) + %epsilon%,0);
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

Solve m2 using lp minimizing w_prowaste;
put log; put 'minimized prowaste mass'/; putclose;
if (m2.modelStat ne 1, abort.noError "no optimal solution found";);
w_prowaste_opt = w_prowaste.l;
* force every pro/ptm waste to be within the optimal limit
v.up(j)$((prowaste(j) or ptmwaste(j)) and S('BIO-protwasted',j) gt 0) = min((w_prowaste_opt / S('BIO-protwasted',j)) + max(%epsilon%,%epsilon%*(w_prowaste_opt / S('BIO-protwasted',j))),v.up(j));
equation m2opt; m2opt.. w_prowaste =l= w_prowaste_opt + max(%epsilon%,w_prowaste_opt*%epsilon%);

model m3 /m2, m2opt/;
m3.optfile=1;
m3.solvelink=5;

Solve m3 using lp minimizing kappEstSlackSum;
put log; put 'minimized kappEstSlackSum'/; putclose;
if (m3.modelStat ne 1, abort.noError "no optimal solution found";);
kappEstSlackSum_opt = kappEstSlackSum.l;
* multipliers just used for tolerance
equation m3opt; m3opt.. kappEstSlackSum =l= kappEstSlackSum_opt + max(%epsilon%, kappEstSlackSum_opt*%epsilon%, card(rxns_with_enz)*%epsilon%)*15;

model m4 /m3, m3opt/;
m4.optfile=1;
m4.solvelink=5;
	
Solve m4 using lp minimizing fluxSlack;
put log; put 'minimized fluxSlack'/; putclose;
if (m4.modelStat ne 1, abort.noError "no optimal solution found";);
fluxSlack_opt = fluxSlack.l;
*fluxSlack.up = fluxSlack_opt * (1+%epsilon%) * card(sm_j_measured);
fluxSlack.up = fluxSlack_opt + (10 * %epsilon% * card(sm_j_measured));
*fluxSlack.up = fluxSlack.l * (1+%epsilon%);
* force each variable in sum to not exceed current sum (+ some leeway for numerical precision) or its prior bounds (to prevent this from adding flexibility where unexpected), to improve solving times and reduce degeneracy
s_v_exp_ub.up(sm_j_measured) = min(fluxSlack.up, s_v_exp_ub.up(sm_j_measured));
s_v_exp_lb.up(sm_j_measured) = min(fluxSlack.up, s_v_exp_lb.up(sm_j_measured));

Solve m4 using lp minimizing fluxSum;
put log; put 'minimized fluxSum'/; putclose;
if (m4.modelStat ne 1, abort.noError "no optimal solution found";);
* force total flux to previous value
fluxSum_opt = fluxSum.l;
fluxSum.up = fluxSum_opt+(10*%vmin%*card(rxns_metab));

ff.nr = 2; put ff;
put m4.modelStat/;
putclose;

file ff2 /%system.FN%.objectives.txt/;
ff2.nr = 2; put ff2; ff2.pc=6;
put 'prosynSlackSum',prosynSlackSum.l:0:11/;
put 'fluxSlackSum',(fluxSlack.l/%nscale%):0:11/;
put 'fluxSum_j_NP',(fluxSum_j_NP.l/%nscale%):0:11/;
put 'fluxSum',(fluxSum.l/%nscale%):0:11/;
putclose;

$include %model_root_path%GAMS/extra_outputs.gms

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

file ffpronodata /%system.FN%.pro_no_data.txt/;
ffpronodata.nr=2; put ffpronodata;
put '/'/;
loop(pro$(not pro_measured(pro)),
	loop(j$pro_prosyn(pro,j),
		if (v.l(j) gt %vmin%,
			put j.tl:0,system.tab,(v.l(j)/%nscale%):0:15/;
		else
			put j.tl:0/;
		);
	);
);
put '/';
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

file ff8 /%system.FN%.enzload_used.txt/;
ff8.nr=2; put ff8;
loop(j$rxns_enzload_used(j),
    put j.tl:0, system.tab, v.l(j):0:15/;
);
putclose;

file ff8a /%system.FN%.kappEstSlack.txt/;
ff8a.nr=2; ff8a.pc=6; put ff8a;
put 'j','v','v_enzload','+ slack','- slack'/;
loop(rxns_with_enz,
	if ( (v.l(rxns_with_enz) gt %vmin% and (kappEstSlackPos.l(rxns_with_enz) gt 0 or kappEstSlackNeg.l(rxns_with_enz) gt 0)),
		put rxns_with_enz.tl:0, v.l(rxns_with_enz):0:15, sum(j1$enzload_rxn_coupling(j1,rxns_with_enz),v.l(j1)):0:15, kappEstSlackPos.l(rxns_with_enz), kappEstSlackNeg.l(rxns_with_enz)/;
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

parameters
focus(c) / #c 1 /
;
variables
objective_fluxes_sum
objective_dummy
;
equations
objective_fluxes_sum_def, focus_def
;
objective_fluxes_sum_def.. objective_fluxes_sum =e= sum(j$objective_fluxes(j), v(j));
focus_def.. objective_dummy =e= sum(c, focus(c) * membrane_SA(c));

Model prot_allo /m4+objective_fluxes_sum_def+focus_def+membrane_SA_eqn/;
*Model prot_allo /all-EnzCap-minPro-ProwasteLim-prot_mass_fraction_max_eqn-membrane_SA_vivo-v_trans_UB_fx-v_trans_LB_fx-v_trans_AA_sum_def/;
prot_allo.optfile = 1;
prot_allo.solvelink=1;
alias(j_ea,j);
alias(c_ea, c);
* Maximize mass fraction of protein in each mitochondrial compartment, to estimate max possible mitochondrial proteome allocation
* If desired, could only check soluble proteins (i.e., (ord(j_ea.tl,12) ne 109)), since surface area constraints should be more effective for membrane-bound proteins; but calculating all values may help for testing
$onMultiR
parameters prot_mass_fraction_max(j) /#j 0/, membrane_SA_exp(c) /#c 0/;
$offMulti
membrane_SA_exp(c) = max(sum(j$(SA(j) and j_c(j,c)),v.l(j)),membrane_SA_exp(c));
prot_mass_fraction_max(j)$j_bio_prot(j) = max(v.l(j),prot_mass_fraction_max(j));
loop(j_ea$(j_bio_prot(j_ea)),
	objective_fluxes(j) = no;
	objective_fluxes(j_ea) = yes;
	Solve prot_allo using lp maximizing objective_fluxes_sum;
	if (prot_allo.modelStat ne 1,
		prot_allo.solvelink=1;
		Solve prot_allo using lp maximizing objective_fluxes_sum;
	);
	if (prot_allo.modelStat eq 1,
		membrane_SA_exp(c) = max(membrane_SA.l(c),membrane_SA_exp(c));
		prot_mass_fraction_max(j)$j_bio_prot(j) = max(v.l(j),prot_mass_fraction_max(j));
	else
		put log; put "*modelstat:" prot_allo.modelStat:0 " *'" j_ea.tl:0 "' " (objective_fluxes_sum.l):0:15/; putclose;
	);
	objective_fluxes(j_ea) = no;
);
* Maximize surface area of mitochondrial membrane proteins
loop(c_ea,
	focus(c) = 0;
	focus(c_ea) = 1;
    Solve prot_allo using lp maximizing objective_dummy;
    if (prot_allo.modelStat ne 1,
		prot_allo.solvelink=1;
		Solve prot_allo using lp maximizing objective_dummy;
	);
    if (prot_allo.modelStat eq 1,
		membrane_SA_exp(c) = max(membrane_SA.l(c),membrane_SA_exp(c));
		prot_mass_fraction_max(j)$j_bio_prot(j) = max(v.l(j),prot_mass_fraction_max(j));
	else
		put log; put "*modelstat:" prot_allo.modelStat:0 " *'" c_ea.tl:0 "' " (objective_dummy.l):0:15/; putclose;
    );
	focus(c_ea) = 0;
);

file f_memSA /membrane_SA_exp.txt/;
f_memSA.nr=2; f_memSA.nz=1e-30; put f_memSA;
put '/'/;
loop(c_ea,
	put "'" c_ea.tl:0 "' " (membrane_SA_exp(c_ea)/%nscale%):0:15/;
);
put '/';
putclose;

file f_protmass /prot_mass_fraction_max.txt/;
f_protmass.nr=2; f_protmass.nz=1e-30; put f_protmass;
put '/'/;
loop(j_ea$j_bio_prot(j_ea),
	put "*fraction=" (prot_mass_fraction_max(j_ea)/v.l('BIOSYN-PROTTOBIO')):0:15/ "'" j_ea.tl:0 "' " (prot_mass_fraction_max(j_ea)):0:15/;
);
put '/';
putclose;
