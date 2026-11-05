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
problematic_rxns(j) /
$if exist problematic_rxns.txt $include problematic_rxns.txt
/
enzload_rxn_coupling_temp(j,j)
;

enzload_rxn_coupling_temp(j,j1)$enzload_rxn_coupling(j,j1) = yes;
*loop(j$problematic_rxns(j),
*	loop(j1$enzload_rxn_coupling_temp(j1,j),
*		enzload_rxn_coupling_temp(j1,j) = no;
*	);
*);
rxns_enzload_nonused(j) = rxns_enzload(j) and not rxns_enzload_used(j);

Parameters
kapp(j)
$include "./kapp_init.txt"
prosynSlackSum_opt
fluxSlack_opt
fluxSum_j_NP_opt
fluxSum_opt
kappEstSlackSum_opt
w_prowaste_opt
;

binary variable y(j);

Variables
fluxSum_j_NP_nonEss
fluxSum_enzload_nonused
y_sum
;

y.fx(j) = 0;
y.up(problematic_rxns) = 1;
y_sum.up = card(problematic_rxns) * 2;

** Disable enzyme load network for reactions that can't be used
*v.fx(j)$rxns_enzload(j) = 0;
* force production of all enzymes involved in used rxns
*v.lo(j)$rxns_enzload_used(j) = %vmin% * (1+1e-8); 
v.up(j)$rxns_enzload_used(j) = %vmax% * %nscale%;

Equations
flux_j_NP_NonEss, flux_enzload_nonused, kappEst(j,rxns_with_enz), kappMaxEst(j,rxns_with_enz), y_sum_def, min_flux_EA(j) 
;

flux_j_NP_NonEss..				fluxSum_j_NP_nonEss =e= sum(j$rxns_NP_nonessential(j), v(j));
flux_enzload_nonused..		fluxSum_enzload_nonused =e= sum(j$rxns_enzload_nonused(j), v(j));
kappEst(j,rxns_with_enz)$enzload_rxn_coupling_temp(j,rxns_with_enz).. (kappEstSlackPos(rxns_with_enz) - kappEstSlackNeg(rxns_with_enz) + v(rxns_with_enz)) * mu =e= v(j) * kapp(rxns_with_enz);
*kappMaxEst(j,rxns_with_enz)$enzload_rxn_coupling_temp(j,rxns_with_enz).. v(rxns_with_enz) * mu =l= v(j) * %kapp_max%;
kappMaxEst(j,rxns_with_enz)$enzload_rxn_coupling_temp(j,rxns_with_enz).. v(rxns_with_enz) * mu =l= (v(j) * %kapp_max%) + (y(rxns_with_enz) * %kapp_max% * 1e99);
y_sum_def.. y_sum =e= sum(j,y(j));
min_flux_EA(j)$rxns_enzload_used(j).. v(j) =g= %vmin% * (1-sum(j1$enzload_rxn_coupling(j,j1),y(j1)));
$include "./enz_alloc_constraints.txt"

* minimize disagreement with proteomics data, while allowing some where needed (e.g., measurement errors)
* exclude v_trans equations since they're not needed and may slow down run times; add back if using them here
Model kapp_calc /all-EnzCap-minPro-ProwasteLim-prot_mass_fraction_max_eqn-membrane_SA_vivo-membrane_SA_eqn-v_trans_UB_fx-v_trans_LB_fx-v_trans_AA_sum_def/;
kapp_calc.optfile = 1;
*Solve kapp_calc using mip minimizing prosynSlackSum;
put log; put 'minimized prosynSlackSum'/; putclose;

* initially allow some minor slack (to account for measurement error) on just proteins w/o modeled functions, to minimize the potential impact on enzymes
* do positive 1st, to minimize impact on ribosome calculations
if (kapp_calc.modelStat ne 1, 
    put log; put 'Model infeasible - relaxing proteomics constraints'/; putclose;
	prosynSlackUB.up(pro_measured)$(pro_wmf(pro_measured) and not unlimited_slack(pro_measured)) = 0.05;
*    Solve kapp_calc using mip minimizing prosynSlackSum;
);
* Try max possible positive slack for proteins w/o modeled functions, to avoid excluding proteins from consideration as ribo subunits
if (kapp_calc.modelStat ne 1, 
	prosynSlackUB.up(pro_measured)$pro_wmf(pro_measured) = inf;
*	Solve kapp_calc using mip minimizing prosynSlackSum;
);
* Try adding small slack to all proteins, to account for measurement error
if (kapp_calc.modelStat ne 1, 
	prosynSlackUB.up(pro_measured)$((not pro_wmf(pro_measured)) and (not unlimited_slack(pro_measured))) = 0.05;
	put log; put 'Allowing small slack for all proteins'/; putclose;
*	Solve kapp_calc using mip minimizing prosynSlackSum;
	put log; put 'Modelstat='kapp_calc.modelstat/; putclose;
);
* Allow max positive slack to all proteins
if (kapp_calc.modelStat ne 1, 
	prosynSlackUB.up(pro_measured) = inf;
	put log; put 'Allowing slack for all proteins'/; putclose;
	Solve kapp_calc using mip minimizing prosynSlackSum;
	put log; put 'Modelstat='kapp_calc.modelstat/; putclose;
);
* Allow negative slack
if (kapp_calc.modelStat ne 1, 
* set bounds to min. height needed to cover full range w/o being excessive (to reduce solving time)
	prosynSlackLB.up(pro_measured) = 1 + %bi_tol%;
*	Solve kapp_calc using mip minimizing prosynSlackSum;
	put log; put 'Modelstat='kapp_calc.modelstat/; putclose;
);

bi_UB=y_sum.up;
file ff_y /y_j.txt/;
while(bi_UB - bi_LB > 1,
	put log; put 'Bisection: ' bi_LB ' to ' bi_UB/; putclose;
	Solve kapp_calc using mip minimizing prosynSlackSum;
	put log; put 'Modelstat='kapp_calc.modelstat/; putclose;
	if ((kapp_calc.modelstat eq 1) or (kapp_calc.modelstat eq 2) or (kapp_calc.modelstat eq 7) or (kapp_calc.modelstat eq 8),
		put ff_y;
		loop(j$(y.l(j) eq 1),
			put j.tl:0/;
		);
		putclose;
		bi_UB=y_sum.l;
	else
		bi_LB = y_sum.up;
	);
	y_sum.up = (bi_UB + bi_LB) / 2;
);

*Solve kapp_calc using mip minimizing y_sum;

if (kapp_calc.modelStat ne 1, put ff; put kapp_calc.modelStat/; putclose; abort.noError "no optimal solution found even after relaxing constraints";);

prosynSlackSum.up = prosynSlackSum.l * (1+%epsilon%);

Solve kapp_calc using mip minimizing w_prowaste;
put log; put 'minimized prowaste mass'/; putclose;
if (kapp_calc.modelStat ne 1, abort.noError "no optimal solution found";);
w_prowaste.up = w_prowaste.l * (1+%epsilon%);

Solve kapp_calc using mip minimizing kappEstSlackSum;
put log; put 'minimized kappEstSlackSum'/; putclose;
if (kapp_calc.modelStat ne 1, abort.noError "no optimal solution found";);
kappEstSlackSum.up = kappEstSlackSum.l * (1+%epsilon%);

Solve kapp_calc using mip minimizing fluxSlack;
put log; put 'minimized fluxSlack'/; putclose;
if (kapp_calc.modelStat ne 1, abort.noError "no optimal solution found";);
fluxSlack.up = fluxSlack.l * (1+%epsilon%);

Solve kapp_calc using mip minimizing fluxSum_j_NP;
put log; put 'minimized fluxSum_j_NP'/; putclose;
if (kapp_calc.modelStat ne 1, abort.noError "no optimal solution found";);
* force rxns that were turned off to stay off
*v.fx(j)$(rxns_with_no_prodata(j) and (v.l(j) eq 0)) = 0;

fluxSum_j_NP.up = fluxSum_j_NP.l * (1+%epsilon%);

*loop(enzload_rxn_coupling(j1,j),
*	if (v.l(j) eq 0, 
*			v.fx(j1) = 0;
*		);
*);

Solve kapp_calc using mip minimizing fluxSum;
put log; put 'minimized fluxSum'/; putclose;
if (kapp_calc.modelStat ne 1, abort.noError "no optimal solution found";);
* force total flux to previous value
fluxSum.up = fluxSum.l*(1+%tol%);

* Solve again, minimizing usage of ENZLOAD fluxes for rxns that aren't being used, to avoid arbitrary reductions in ENZLOAD fluxes that can increase kapps by reducing the denominator; how to fix this is unclear.
Solve kapp_calc using mip minimizing fluxSum_enzload_nonused;
put log; put 'minimized flux through ENZLOAD rxns that aren''t being used'/; putclose;
if (kapp_calc.modelStat ne 1, abort.noError "no optimal solution found";);
fluxSum_enzload_nonused.up = fluxSum_enzload_nonused.l * (1+%epsilon%);

* Solve again, encouraging more equal use of all enzymes
* NOTE: disabled this step since it can lead to arbitrary reductions in ENZLOAD fluxes, even when other ones aren't being used. This can increase kapps by reducing the denominator; how to fix this is unclear.
*Solve kapp_calc using mip minimizing slackSum;
*put log; put 'minimized uneven enzload distribution'/; putclose;

* Solve again, minimizing all fluxes (not just metabolic ones)
Solve kapp_calc using mip minimizing v_sum;
put log; put 'minimized v_sum'/; putclose;
if (kapp_calc.modelStat ne 1, abort.noError "no optimal solution found";);
*v_sum.up = v_sum.l * (1+%epsilon%);

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


file ff8 /%system.FN%.enzload_used.txt/;
ff8.nr=2; put ff8;
loop(j$rxns_enzload_used(j),
    put j.tl:0, system.tab, v.l(j):0:15/;
);
putclose;

file ff8a /%system.FN%.kappEstSlack.txt/;
ff8a.nr=2; ff8a.pc=6; put ff8a;
put 'j','v','v_enzload','+ slack','- slack'/;
loop(enzload_rxn_coupling(j1,rxns_with_enz),
	if ( (v.l(rxns_with_enz) gt %vmin% and (kappEstSlackPos.l(rxns_with_enz) gt 0 or kappEstSlackNeg.l(rxns_with_enz) gt 0)),
		put rxns_with_enz.tl:0, v.l(rxns_with_enz):0:15, v.l(j1):0:15, kappEstSlackPos.l(rxns_with_enz), kappEstSlackNeg.l(rxns_with_enz)/;
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

Model prot_allo /all-EnzCap-minPro-ProwasteLim-prot_mass_fraction_max_eqn-membrane_SA_vivo/;
* Maximize mass fraction of protein in each mitochondrial compartment, to estimate max possible mitochondrial proteome allocation
* If desired, could only check soluble proteins (i.e., (ord(j_ea.tl,12) ne 109)), since surface area constraints should be more effective for membrane-bound proteins; but calculating all values may help for testing
file f_protmass /prot_mass_fraction_max.txt/;
f_protmass.nr=2; f_protmass.nz=1e-30; put f_protmass;
alias(j_ea,j);
alias(c_ea, c);
put '/'/;
loop(j_ea$(j_bio_prot(j_ea)),
	loop(c_ea$(j_c(j_ea,c_ea)),
        objective_fluxes(j) = no;
        objective_fluxes(j_ea) = yes;
		Solve prot_allo using mip maximizing objective_fluxes_sum;
		put "'" j_ea.tl:0 "' " (objective_fluxes_sum.l/v.l('BIOSYN-PROTTOBIO')):0:15/;
        objective_fluxes(j_ea) = no;
    );
);
put '/'/;
putclose;

* Maximize surface area of mitochondrial membrane proteins
file f_memSA /membrane_SA_exp.txt/;
f_memSA.nr=2; f_memSA.nz=1e-30; put f_memSA;
put '/'/;
loop(c_ea,
	focus(c) = 0;
	focus(c_ea) = 1;
    Solve prot_allo using mip maximizing objective_dummy;
    put "'" c_ea.tl:0 "' " membrane_SA.l(c_ea):0:15/;
	focus(c_ea) = 0;
);
put '/'/;
putclose;
