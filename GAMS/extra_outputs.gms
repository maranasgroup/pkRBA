* extra output files used after solving a model

file fi_v_sm / %gms%.SM_flux.txt /;
put fi_v_sm;
loop(sm_j,
	put sm_j.tl:0, system.tab, (sum(j$dir(sm_j,j), dir(sm_j,j)*v.l(j)) / %nscale%):0:15, system.tab, (v_sm_lo(sm_j)/%nscale%), system.tab, (v_sm_up(sm_j)/%nscale%)/;
);

file fi_v_bounds / %gms%.flux_with_bounds.txt /;
put fi_v_bounds;
put 'j',system.tab,'flux',system.tab,'LB',system.tab,'UB'/;
loop(j,
	put j.tl:0, system.tab, (v.l(j) / %nscale%):0:15, system.tab, (v.lo(j)/%nscale%):0:15, system.tab, (v.up(j)/%nscale%):0:15/;
);
putclose;
