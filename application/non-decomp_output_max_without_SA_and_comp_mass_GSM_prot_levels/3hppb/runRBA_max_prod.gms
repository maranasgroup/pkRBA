* run RBA without bisection, using default settings and outputs
$if not setGlobal gms $setGlobal gms %system.FN%
* turns off stopping caused by domain violations; useful for testing decompartmentalized and non-decompartmentalized models in 1 place, but best left off in most other cases as quality control 
$onWarning

$include "./runRBA_GAMS_settings.txt"

* if initial solve is infeasible, try adjustments (e.g., changing growth rate)
$setGlobal adjust_constraints_if_infeas 0 
* % of min experimental yield model must at least match
$setGlobal percent_min_exp_yield 0
* default carbon slack (disabled by default)
$setGlobal carbonSlack 0

$include %model_root_path%GAMS/defaults.gms

$setGlobal min_media_path %null_path%
$set path "%gms_path%RBA_rxns_EXREV_minimal.txt"
$if exist %path% $setGlobal min_media_path %path%
sets
rxns_add(j) "list of heterologous rxns added for model application"
$include "%rxns_add_path%"
min_media(j)
$include %min_media_path%
undesired(j) "rxns to minimize in subsequent solves (such as secretions that may adversly impact titer)"
$ifThen exist %gms_path%decompartmentalized_rxns.txt
$include %gms_path%decompartmentalized_rxns.txt
$else
/ #secretion /
$endIf
;
* remove product from undesired list if it was included, to avoid minimization of the product flux in subsequent solves
undesired('%vprod%') = no;

parameters
mw(j)
$include "yield_denominators.txt"
mw_other_possible_substrates(j)
$include "yield_denominators_expanded.txt"
mu_current
;

* add substrate to undesired list, to encourage higher yields by minimizing substrate use
undesired(j)$mw(j) = yes;

variables
v_prod
mv_sub_sum "sum of (substrate flux * MW) for all substrates"
undesired_sum "sum of rxns to minimize in subsequent solves"
;

* Turn off all heterologous rxns rxns_add
* Corresponding rxns for specific product will be turned on in phenotype.txt file
v.up(j)$rxns_add(j) = 0;

$onMultiR
$include phenotype.txt

equations
mvSubSum, minYield, maxProd, def_undesired_sum
;
mvSubSum.. mv_sub_sum =e= sum(j$mw(j),v(j)*mw(j));
minYield.. minProductYield * mv_sub_sum * %percent_min_exp_yield% / 100 =l= %prod_mw% * v('%vprod%');
maxProd..			v_prod =e= -v('%vprod%');
def_undesired_sum..	undesired_sum =e= sum(undesired(j),v(j));

Model rba
/all/;
rba.optfile = 1;
* save copy of all variable values for later use, in case future solves are infeas. but we want to reference past results anyway
rba.savepoint=2;
$offMulti

*turn this on only when needed to export parameters
*execute_unload 'allparams.gdx';
solve rba using lp minimizing v_prod;

$if not setGlobal yield_abridged $setGlobal yield_abridged 0
* set to 2 to skip extra steps; useful for faster runs where FVA doesn't matter as much and risks causing infeasibilities
$setGlobal solve_once 2
$ifThenE.abridgedtwo (%yield_abridged%+%solve_once%<2)
scalar last_feas /1/;
if((rba.modelstat = 1 or rba.modelstat = 2 or rba.modelstat = 5 or rba.modelstat = 8),
	v_prod.up = v_prod.l + %bi_tol%;
	solve rba using lp minimizing undesired_sum;
* load backup if needed to ensure some values are present
	if(not (rba.modelstat = 1 or rba.modelstat = 2 or rba.modelstat = 5 or rba.modelstat = 8),
		put log; put 'reading past solution'; putclose;
		execute_loadpoint 'rba_p1.gdx';
		rba.solvelink = 1;
	else
		last_feas=2;
		undesired_sum.up = undesired_sum.l*(1+%bi_tol%);
	);
* parsimony
	solve rba using lp minimizing v_sum;
	if(not (rba.modelstat = 1 or rba.modelstat = 2 or rba.modelstat = 5 or rba.modelstat = 8),
		put log; put 'reading past solution'; putclose;
		rba.solvelink = 1;
		if(last_feas=1,
			execute_loadpoint 'rba_p1.gdx'; 
		else
			execute_loadpoint 'rba_p2.gdx'; 
		);
* solve again, minimizing z_dummy to avoid skewing the results differently than the original solve while quickly correcting the modelstat
		solve rba using lp minimizing z_dummy;
	);
	v_sum.up = v_sum.l*(1+%bi_tol%);
);
$endIf.abridgedtwo

* solve again, minimizing mv_sub_sum to find a better yield by discouraging wasteful substrate use 
* NOTE: (get rid of this if linear fractional programming version is developed)
* save current uptake bounds as parameters

$ifThen.abridged %yield_abridged% == 0
Parameters
v_up(j), vprod_max
;
v_up(j) = v.up(j);
vprod_max = v.l('%vprod%');
v.lo('%vprod%') = vprod_max;
*if((rba.modelstat eq 1),
*	Solve rba using lp minimizing mv_sub_sum;
*);

* If not optimal, try again with bounds on vprod relaxed in case it's demanding too much.
* This isn't done initially in case setting specific bounds is important for replicating experimental results (e.g., for making multiple products).
v.up('%vprod%') = 1e3 * %nscale%;
v.lo('%vprod%') = 0 * %nscale%;

* If not optimal, find lowest uptakes that still allow optimal growth and exp. yields
* Retry while relaxing uptake upper bounds via slack variables, then minimize their sum
* Included in case uptake rate estimates are inaccurate.

* allow all possible uptakes in medium to reach max flux
v.up(j)$(uptake(j) and v_up(j) gt 0) = 1e3 * %nscale%;

Variables 
slack(j), slacksum
;
* default
slack.fx(j) = 0;
* allow slack for uptakes
slack.up(j)$(uptake(j) and v.up(j) gt 0) = inf;
slack.lo(j)$(uptake(j) and v.up(j) gt 0) = 0 * %nscale%;
* use initial uptake bounds as slack for other reactions
*$include "RBA_GAMS_defaults_from_FBA_initial.txt"
Equations slackcons, slackobj;
slackcons(j)$(slack.up(j) ne 0).. v(j) =l= v_up(j) + slack(j);
slackobj.. slacksum =e= sum(j, slack(j));
Model rba_slack /all/;

* added in case infeasibilities occur, and FBA predictions should be adjusted to match
$onMulti
Model fba /
Obj, Stoic, mvSubSum, minYield, slackcons, slackobj
/;

while((rba.modelstat ne 1 and (%adjust_constraints_if_infeas% ne 0)),
	rba_slack.optfile = 1;
	display 'UPDATE: find lowest uptakes that still allow optimal growth and exp. yields';
	Solve rba_slack using lp minimizing slacksum;
	if(rba_slack.modelstat eq 1, break;);

* If still not optimal, suggests minYield is wrong or the model's parameters don't perfectly match exp. conditions
* Try again, but find the highest yields under the original conditions
* turn off yield requirements
	minProductYield = 0;
	display 'UPDATE: find highest yields under the original conditions';
	Solve rba_slack using lp minimizing z;
	if(rba_slack.modelstat eq 1, break;);

* if still infeas., test if any media composition would work
	v.up(j)$(uptake(j)) = 1e3 * %nscale%;
	slack.up(j)$(uptake(j)) = 1e10 * %nscale%;
* find smallest uptake slacks needed for near-maximum flux
	display 'UPDATE: find smallest uptakes needed for near-maximum product flux';
	Solve rba_slack using lp minimizing slacksum;
	if(rba_slack.modelstat eq 1,
		v.fx("%vprod%") = v.l("%vprod%") * 0.99;
		Solve rba_slack using lp minimizing z;
	);
	break
);

$include %model_root_path%GAMS/default_outputs.gms
$include %model_root_path%GAMS/extra_outputs.gms

file f_undes /runRBA.undesired_rxns_used.txt/;
f_undes.pc=6;
put f_undes;
put 'rxn','flux','LB','UB','v_sm_lo','v_sm_up'/;
loop(j$(undesired(j) and v.l(j) ne 0),
	put j.tl:0, (v.l(j)/%nscale%):0:15,(v.lo(j)/%nscale%):0,(v.up(j)/%nscale%):0,;
	loop(sm_j$dir(sm_j,j),
		put (v_sm_lo(sm_j)/%nscale%):0,(v_sm_up(sm_j)/%nscale%):0;
		break;
	);
	put /;
);

* write slacks and uptakes to file
file ff3 /runRBA.slack.txt/;
ff3.pc=6;
put ff3;
put 'rxn','flux','slack'/;
loop(j$(uptake(j) and v.up(j) gt 0),
	put j.tl:0, (v.l(j) / %nscale%):0:15,slack.l(j):0:15/;
);
putclose ff3;

if(rba.modelstat eq 1,
	abort.noError 'Optimal solution found';
);
if((%adjust_constraints_if_infeas% eq 0),
	abort.noError 'no optimal solution; set adjust_constraints_if_infeas to 1 to try tweaks'
);
* find smallest uptakes on minimal media
* set all other uptakes to 0
v.up(j)$uptake(j) = 0;
v.up(j)$(uptake(j) and min_media(j)) = 1e3 * %nscale%;
display 'UPDATE: find smallest uptakes needed for min exp. growth+yields on minimal media';
Solve rba_slack using lp minimizing slacksum;

v.lo('%vprod%') = 0;
v.up('%vprod%') = 1e3 * %nscale%;
* force slacks to not exceed RBA levels
slack.up(j) = slack.l(j);
* force uptakes to not exceed RBA levels
v.up(j)$uptake(j) = v.l(j);
display 'UPDATE: find FBA fluxes on minimal media using prior uptakes';
solve fba using lp minimizing z;
file ffba /runRBA.FBA_rerun_fluxes.txt/;
put ffba;
ffba.pc=6;
loop(j,
	if ( (v.l(j) gt 0),
		put j.tl:0, 'v', (v.l(j) / %nscale%):0:15/;
	);
);
putclose ffba;

* output yield

set n "percent of growth rate from prior solutions" /10, 0/;
parameters
yield_denominator
;
file yields /runRBA.yields.txt/;
put yields;
yields.pc=6;
put 'prod', 'method', 'growth rate', 'yield'/;
yield_denominator = sum(j$mw(j),v.l(j)*mw(j));
* if original substrates are unused, assume glucose was used
if(yield_denominator eq 0,
	yield_denominator = sum(j$mw_other_possible_substrates(j),v.l(j)*mw_other_possible_substrates(j));
);
put '%vprod%','fba',(v.l('%bio%') / %nscale%):0:15,(%prod_mw% * v.l('%vprod%') / yield_denominator):0:15/;
* find RBA yield under these conditions
solve rba using lp minimizing z;
yield_denominator = sum(j$mw(j),v.l(j)*mw(j));
* if original substrates are unused, assume glucose was used
if(yield_denominator eq 0,
	yield_denominator = sum(j$mw_other_possible_substrates(j),v.l(j)*mw_other_possible_substrates(j));
);
put '%vprod%','rba',(v.l('%bio%') / %nscale%):0:15,(%prod_mw% * v.l('%vprod%') / yield_denominator):0:15/;
* allow PROWASTE rxns since without growth, there's no protein sink
v.up(j)$prowaste(j) = 1e3 * %nscale%;

loop(n,
	mu_current = mu * n.val/100;
	v.fx('%bio%') = mu_current * %nscale%;
	solve fba using lp minimizing z;
	yield_denominator = sum(j$mw(j),v.l(j)*mw(j));
* if original substrates are unused, assume glucose was used
	if(yield_denominator eq 0,
		yield_denominator = sum(j$mw_other_possible_substrates(j),v.l(j)*mw_other_possible_substrates(j));
	);
	put '%vprod%','fba',(v.l('%bio%') / %nscale%):0:15,(%prod_mw% * v.l('%vprod%') / yield_denominator):0:15/;
	solve rba using lp minimizing z;
	yield_denominator = sum(j$mw(j),v.l(j)*mw(j));
* if original substrates are unused, assume glucose was used
	if(yield_denominator eq 0,
		yield_denominator = sum(j$mw_other_possible_substrates(j),v.l(j)*mw_other_possible_substrates(j));
	);
	put '%vprod%','rba',(v.l('%bio%') / %nscale%):0:15,(%prod_mw% * v.l('%vprod%') / yield_denominator):0:15/;
);

putclose yields;

file frba2 /runRBA.RBA_last_rerun_fluxes.txt/;
put frba2;
frba2.pc=6;
loop(j,
	if ( (v.l(j) gt 0),
		put j.tl:0, 'v', (v.l(j) / %nscale%):0:15/;
	);
);
putclose frba2;
$endIf.abridged

