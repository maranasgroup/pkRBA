* Run FVA (Flux Variability Analysis) on the RBA model
* To use in conjunction with another script, set the objective function to its optimal value from that script and then include this file
$if not setGlobal gms $setGlobal gms %system.FN%

* turn off unnecessary rxns to improve solving time
$ifThen %gms%==FVA
$INLINECOM /*  */
$include ../paths.txt
$offInline
* run RBA without bisection, using default settings and outputs
$include %model_root_path%GAMS/defaults.gms
$endIf


Sets
rxns_to_test(j) "rxns to include in FVA" /#j
/
dummy(j) "dummy set for iterating over rxns"
;

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

rxns_to_test(j)$((v.lo(j) gt -%vmin% and v.up(j) lt %vmin%) or (v.lo(j) eq v.up(j))) = no;
alias(j2,j);

loop(j2$(rxns_to_test(j2)),
    dummy(j)$(ord(j) ne ord(j2)) = no;
    dummy(j)$(ord(j) eq ord(j2)) = yes;

	high(j2) = v.lo(j2) / %nscale%;
	low(j2) = v.lo(j2) / %nscale%;
	high_stat(j2) = -1;
	low_stat(j2) = -1;
	Solve rba using lp maximizing z_fva;
	high_stat(j2) = rba.modelstat;
	If((rba.modelstat eq 1) and (z_fva.l gt v.lo(j2)),
		high(j2) = z_fva.l / %nscale%;
		solve rba using lp minimizing z_fva;
		low_stat(j2) = rba.modelstat;
		if(rba.modelstat eq 1,
			low(j2) = z_fva.l / %nscale%;
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
v.up(j)$(high_stat(j) eq 1) = high(j) * %nscale%;
v.lo(j)$(low_stat(j) eq 1) = low(j) * %nscale%;
loop(j2,
	put j2.tl:66, v.lo(j2):15:8, v.up(j2):15:8/;
);
putclose f_fva3;

$include %model_root_path%GAMS/default_outputs.gms