$if not setGlobal gms $setGlobal gms %system.FN%
$include enz_alloc_setup.gms
put log;
loop(j,
	put j.tl:0 ' '  v.lo(j):0 ' ' v.up(j):0/;
);
putclose;
if(smax(pro_measured, prosynSlackUB.up(pro_measured) + prosynSlackLB.up(pro_measured)) le %vmin%,
	loop(pro_measured,
		if(sum(j$(pro_prosyn(pro_measured,j) and v.up(j) gt 0),1) eq 1,
			v.fx(j)$(pro_prosyn(pro_measured,j) and v.up(j) gt 0) = v.l(j);
		);
	);
);
$include FVA_turbo.gms
