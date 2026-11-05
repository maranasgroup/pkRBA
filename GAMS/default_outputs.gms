* default output files used after solving a model; see defaults.gms for name details
put ff;
put rba.modelStat/;
putclose;

put fi_v;
*TODO: experiment w/ using scientific notation for numbers to ensure that very small numbers are not rounded to 0 in the output files
loop(j,
	if ( (v.l(j) gt 0),
		put j.tl:0, system.tab, 'v', system.tab, (v.l(j) / %nscale%):0:15/;
	);
);
putclose;

put fi_vgs;
*TODO: experiment w/ using scientific notation for numbers to ensure that very small numbers are not rounded to 0 in the output files
loop(j,
	if ( (v.l(j) gt 0),
		put j.tl:0, system.tab, 'v', system.tab, (v.l(j)):0:15/;
	);
);
putclose;
