$include enz_alloc_setup.gms
put log;
loop(j,
	put j.tl:0 ' '  v.lo(j):0 ' ' v.up(j):0/;
);
putclose;
$include FVA.gms
