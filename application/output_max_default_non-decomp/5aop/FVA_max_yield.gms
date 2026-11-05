* find max yields, do FVA
$setGlobal gms %system.FN%
$setGlobal yield_abridged 0
* set to 1 to test each flux extreme, 0 to not test it
*TODO: Get test_max working someday
scalar test_min /0/
test_max /1/;
$include runRBA_max_prod.gms
$onMulti

Sets
rxns_tested(j) "rxns to include in FVA"
dummy(j) "dummy set for iterating over rxns"
;
$offMulti
if ((card(rxns_tested) = 0),
	rxns_tested(j) = yes;
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
model rba /all/;
rba.optfile = 1;
* stop saving backups
rba.savepoint=0;
$offMulti
*Updates output while solving, so errors won't get rid of output
file f_fva /%gms%.txt/;
* Clear old file
put f_fva;
put "rxn":66, "min_flux":15, "max_flux":15, "high_stat":10, "low_stat":10/;
* Append old file
f_fva.ap=1;
f_fva.nr=2;

*$if exist %phenotype_path% $include %phenotype_path%

rxns_tested(j)$(v.lo(j) gt -%vmin% and v.up(j) lt %vmin% or (v.lo(j) eq v.up(j)) or (sum(sm_j$dir(sm_j,j),1) ge 1 and sum(sm_j$dir(sm_j,j),v_sm_lo(sm_j)-v_sm_up(sm_j)) eq 0)) = no;
alias(j2,j);

loop(j2$(rxns_tested(j2)),
    dummy(j)$(ord(j) ne ord(j2)) = no;
    dummy(j)$(ord(j) eq ord(j2)) = yes;
    display 'Performing FVA for rxn', j2;

	high(j2) = v.lo(j2) / %nscale%;
	low(j2) = v.lo(j2) / %nscale%;
	high_stat(j2) = -1;
	low_stat(j2) = -1;
	Solve rba using lp maximizing z_fva;
	high_stat(j2) = rba.modelstat;
	If((rba.modelstat eq 1) and (z_fva.l gt v.lo(j2)),
		high(j2) = z_fva.l / %nscale%;
		if(test_min,
			solve rba using lp minimizing z_fva;
			low_stat(j2) = rba.modelstat;
			if(rba.modelstat eq 1,
				low(j2) = z_fva.l / %nscale%;
			);
		);
	);

	put f_fva;
	put j2.tl:66, low(j2):15:8, high(j2):15:8, high_stat(j2):10:0, low_stat(j2):10:0/;
	putclose f_fva;
);

file f_fva2 /%gms%.bounds.txt/;
f_fva2.ap=1;
f_fva2.nr=2;
put f_fva2;
put "":66, "min_flux":15, "max_flux":15/;
loop(j2$(rxns_tested(j2)),
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
v.up(j)$(high_stat(j) eq 1) = high(j) * %nscale%;
v.lo(j)$(low_stat(j) eq 1) = low(j) * %nscale%;
loop(j2,
	put j2.tl:66, v.lo(j2):15:8, v.up(j2):15:8/;
);
putclose f_fva3;

$include %model_root_path%GAMS/default_outputs.gms