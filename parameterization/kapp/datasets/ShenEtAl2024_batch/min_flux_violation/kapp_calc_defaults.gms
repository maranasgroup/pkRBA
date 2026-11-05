* everything prior to model declaration for use in other files (e.g., min_flux_violation.gms, enz_alloc.gms)
* Created to help ensure consistency between kapp calculation steps
$INLINECOM /*  */

* whether to ignore blocked and blocked_SM rxns (no by default)
$if not setGlobal ignore_blocked $setGlobal ignore_blocked 0
$if not setGlobal ignore_blocked_j $setGlobal ignore_blocked_j %ignore_blocked%
$if not setGlobal ignore_blocked_sm $setGlobal ignore_blocked_sm %ignore_blocked%

* ignore constraints requiring production of measured but unused proteins from kapp calculations
$setGlobal ignore_kapps 1
* allow proteins w/o modeled functions to be used; disable to reduce solver time
$setGlobal ignore_pwmf 0

$include ../paths.txt
* slacks turned off by default, 
* 	but included to account for measurement errors when needed
$setGlobal prosynSlackAllow 0
* small value needed to ensure sequential problems aren't infeasible due to rounding errors
$setGlobal epsilon 1e-4
* update to match solver tolerance
$setGlobal tol 1e-2

$include "./min_flux_violation_GAMS_settings.txt"
* important to avoid errors when reading defaults file
$offInline

$include %model_root_path%GAMS/defaults.gms

v.up(j)$prowaste(j) = %vmax% * %nscale%;

Sets
rxns_with_no_prodata(j) "enzymatic rxns without proteomics data for all their enz subunits"
$include "%rxns_with_no_prodata_path%"
rxns_metab(j)
$include "%rxns_metab_path%"
unlimited_slack(pro) "proteins allowed to have unlimited slack."
rxns_enzsyn(j)
$include "%rxns_enzsyn_path%"
rxns_enzload(j)
$include "%rxns_enzload_path%"
enzload_with_no_prodata(j) "enzload for enzymatic rxns without proteomics data for all their enz subunits"
$include "./enzload_enz_with_no_prodata.txt"
sm_j_lumped "rxns representing combinations of stoich. model rxns"
$include "../sm_j_lumped.txt"
pro_wmf(pro) /
$set path %model_root_path%build_model/proteins_without_modeled_functions.txt
$if exist %path% $include %path%
/
pro_measured(pro) "proteins with measured proteomics data"
$include "%pro_measured_path%"
sm_j_measured(sm_j) "stoich. model rxns with measured flux data"
$include "%sm_j_measured_path%"
blocked(j) /
$if exist "blocked_rxns.txt" $include "blocked_rxns.txt"
/
blocked_sm(sm_j)
$ifThen exist "blocked_rxns_SM.txt"
$include "blocked_rxns_SM.txt"
$else $include %null_path%
$endIf 
;
if(%ignore_blocked_j% eq 0,
	v.fx(j)$blocked(j)=0;
);
if(%ignore_blocked_sm% eq 0,
	v_sm_lo(sm_j)$blocked_sm(sm_j)=0;
	v_sm_up(sm_j)$blocked_sm(sm_j)=0;
);
*As in scRBA and rtRBA, allow slack to be used freely for ribosome subunits, to account for missing subunit paralogs or ambiguous stoichiometry.
unlimited_slack(pro)$sum(r,pro_ribo(r,pro)) = yes;

Parameters
v_exp_pro(pro)
$include "%v_exp_pro_path%"
v_exp_lb(sm_j)
$include "%v_exp_lb_path%"
v_exp_ub(sm_j)
$include "%v_exp_ub_path%"
sm_j_lumped_mappings(sm_j_lumped,sm_j) "grouping stoich. model rxns to their respective lumped counterparts"
$include "../sm_j_lumped_mappings.txt"
;

* improve solving time by removing proteins w/o modeled functions
if(%ignore_pwmf% eq 1,
	loop(pro$(pro_wmf(pro) and (not pro_measured(pro))),
		loop(j$pro_prosyn(pro,j),
			v.fx(j) = 0;
			loop(j1$pro_syn_waste(j,j1),
				v.fx(j1) = 0;
			);
		);
	);
);

* slacks for allowing fluxes to deviate from measured values when necessary
Variables
w_prowaste, prosynSlackSum, kappEstSlackSum, fluxSum_j_NP, fluxSum, fluxSlack, s_v_exp_lb(sm_j_measured), s_v_exp_ub(sm_j_measured), s_v_exp_lb_lumped(sm_j_lumped), s_v_exp_ub_lumped(sm_j_lumped), prosynSlackLB(pro_measured), prosynSlackUB(pro_measured), EnzLoadSlackPos(rxns_with_enz), EnzLoadSlackNeg(rxns_with_enz), kappEstSlackPos(rxns_with_enz), kappEstSlackNeg(rxns_with_enz), slackSum
;
w_prowaste.lo = 0;
prosynSlackSum.lo = 0;
kappEstSlackSum.lo = 0;
fluxSum_j_NP.lo = 0;
fluxSum.lo = 0;
fluxSlack.lo = 0;
slackSum.lo = 0;

prosynSlackLB.lo(pro_measured) = 0; prosynSlackLB.up(pro_measured) = %prosynSlackAllow%;
prosynSlackUB.lo(pro_measured) = 0; prosynSlackUB.up(pro_measured) = %prosynSlackAllow%;
prosynSlackLB.up(pro_measured)$unlimited_slack(pro_measured) = 1;
prosynSlackUB.up(pro_measured)$unlimited_slack(pro_measured) = inf;
* allow changes in either direction; not much more or less
s_v_exp_ub.lo(sm_j_measured) = 0; s_v_exp_ub.up(sm_j_measured) = ((10*%bi_tol%) + abs(v_sm_up(sm_j_measured) - v_exp_ub(sm_j_measured))) * %nscale%;
s_v_exp_lb.lo(sm_j_measured) = 0; s_v_exp_lb.up(sm_j_measured) = ((10*%bi_tol%) + abs(v_sm_lo(sm_j_measured) - v_exp_lb(sm_j_measured))) * %nscale%;
* slacks for lumped fluxes
s_v_exp_lb_lumped.lo(sm_j_lumped) = 0; s_v_exp_lb_lumped.up(sm_j_lumped) = 2 * %vmax% * %nscale%;
s_v_exp_ub_lumped.lo(sm_j_lumped) = 0; s_v_exp_ub_lumped.up(sm_j_lumped) = 2 * %vmax% * %nscale%;
* slacks for enzyme load, to allow for more equal use of enzymes
EnzLoadSlackPos.lo(rxns_with_enz) = 0; EnzLoadSlackPos.up(rxns_with_enz) = 10 * %vmax% * %nscale%;
EnzLoadSlackNeg.lo(rxns_with_enz) = 0; EnzLoadSlackNeg.up(rxns_with_enz) = 10 * %vmax% * %nscale%;
Equations
Weight_Prowaste, kappEstSlack, ProsynSlacks, Obj1, Obj2, Obj3, sm_LB_exp, sm_UB_exp, fluxSlackBounds, lumped_flux(sm_j_lumped)
;

Weight_Prowaste..               w_prowaste =e= v('PROWASTE-TOTALPROTEIN');
kappEstSlack..		kappEstSlackSum =e= sum(rxns_with_enz, kappEstSlackPos(rxns_with_enz) + kappEstSlackNeg(rxns_with_enz));
ProsynSlacks(pro_measured)$(smax(j$pro_prosyn(pro_measured,j),v.up(j)) gt %vmin%).. sum(j$pro_prosyn(pro_measured,j),v(j)) =e= v_exp_pro(pro_measured) * (1 - prosynSlackLB(pro_measured) + prosynSlackUB(pro_measured)) * %nscale%;
*DefSumPTM..			SumPTM =e= sum(j$ptm(j), v(j));
Obj1..				prosynSlackSum =e= sum(pro_measured$(smax(j$pro_prosyn(pro_measured,j),v.up(j)) gt %vmin%), prosynSlackLB(pro_measured) + prosynSlackUB(pro_measured));
Obj2..				fluxSum_j_NP =e= sum(j$rxns_with_no_prodata(j), v(j));
Obj3..				fluxSum =e= sum(j$rxns_metab(j), v(j));
* sm upper and lower bounds for fluxes (if data available); slacks included in case necessary
sm_LB_exp(sm_j_measured)$v_exp_lb(sm_j_measured).. sum(j$dir(sm_j_measured,j),dir(sm_j_measured,j)*v(j)) =g= (v_exp_lb(sm_j_measured) * %nscale%) - s_v_exp_lb(sm_j_measured);
sm_UB_exp(sm_j_measured)$v_exp_ub(sm_j_measured).. sum(j$dir(sm_j_measured,j),dir(sm_j_measured,j)*v(j)) =l= (v_exp_ub(sm_j_measured) * %nscale%) + s_v_exp_ub(sm_j_measured);
fluxSlackBounds..		fluxSlack =e= sum(sm_j_measured, s_v_exp_lb(sm_j_measured) + s_v_exp_ub(sm_j_measured)) + sum(sm_j_lumped, s_v_exp_lb_lumped(sm_j_lumped) + s_v_exp_ub_lumped(sm_j_lumped));
lumped_flux(sm_j_lumped).. sum(sm_j$sm_j_lumped_mappings(sm_j_lumped,sm_j),sum(j$dir(sm_j,j),dir(sm_j,j)*v(j))) =e= (smax(sm_j,sm_j_lumped_mappings(sm_j_lumped,sm_j)) * %nscale%) + s_v_exp_ub_lumped(sm_j_lumped) - s_v_exp_lb_lumped(sm_j_lumped);

*Equation Obj5; Obj5.. slackSum =e= sum(rxns_with_enz, EnzLoadSlackNeg(rxns_with_enz) + EnzLoadSlackPos(rxns_with_enz));

* allows command-line to output text; useful for debugging
file log /''/; 

* check if any pro_measured proteins are blocked; report which ones
put log;
loop(pro_measured$(v_exp_pro(pro_measured) gt 0 and smax(j$pro_prosyn(pro_measured,j),v.up(j)) le %vmin%),
	put 'pro_measured ' pro_measured.tl:0 ' is blocked and cannot be produced, but has measured proteomics data of ' v_exp_pro(pro_measured):0:15/;
);
putclose;
