* Run FVA (Flux Variability Analysis) on the RBA model (over sm_j rxns only)
* To use in conjunction with another script, set the objective function to its optimal value from that script and then include this file
$if not setGlobal gms $setGlobal gms %system.FN%

$INLINECOM /*  */
$include ../paths.txt
$offInline

$include %model_root_path%GAMS/defaults.gms

Sets
rxns_tested(sm_j) "rxns to include in FVA" / #sm_j /
dummy(sm_j) "dummy set for iterating over rxns"
;

Parameters
high(sm_j), low(sm_j), high_stat(sm_j), low_stat(sm_j)
;

Variables
z_fva
;

Equations
Obj_FVA
;

Obj_FVA..           z_fva =e= sum(sm_j$dummy(sm_j),sum(j$dir(sm_j,j),dir(sm_j,j)*v(j)));
$onMultiR
model rba /all/;
rba.optfile = 1;
$offMulti
*Updates output while solving, so errors won't get rid of output
file f_fva /%system.FN%.txt/;
put f_fva;
* Clear old file
put "rxn",system.tab,"min_flux",system.tab,"max_flux",system.tab,"min_flux_modelstat",system.tab,"max_flux_modelstat"/;
* Append old file
f_fva.ap=1;

$if exist %phenotype_path% $include %phenotype_path%

alias(sm_j1,sm_j);
loop(sm_j1$(rxns_tested(sm_j1)),
    dummy(sm_j)$(ord(sm_j) ne ord(sm_j1)) = no;
    dummy(sm_j)$(ord(sm_j) eq ord(sm_j1)) = yes;
    display 'Performing FVA for rxn', sm_j1;

    high(sm_j1) = v_sm_lo(sm_j1);
    low(sm_j1) = v_sm_lo(sm_j1);
    high_stat(sm_j1) = -1;
    low_stat(sm_j1) = -1;
    
    if((v_sm_lo(sm_j1) ne 0) or (v_sm_up(sm_j1) ne 0),
		Solve rba using lp maximizing z_fva;
		high_stat(sm_j1) = rba.modelstat;
        If((rba.modelstat eq 1) and (z_fva.l gt v_sm_lo(sm_j1)),
			high(sm_j1) = z_fva.l / %nscale%;
            solve rba using lp minimizing z_fva;
            low_stat(sm_j1) = rba.modelstat;
            if(rba.modelstat eq 1,
                low(sm_j1) = z_fva.l / %nscale%;
            );
        );
	);
    put f_fva;
    put sm_j1.tl:0:50,system.tab,low(sm_j1):0:8,system.tab,high(sm_j1):0:8,system.tab,low_stat(sm_j1):0:0,system.tab,high_stat(sm_j1):0:0/;
    putclose f_fva
);

$include %model_root_path%GAMS/default_outputs.gms