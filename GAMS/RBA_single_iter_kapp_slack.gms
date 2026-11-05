* run RBA without bisection, using default settings and outputs
$if not setGlobal gms $setGlobal gms %system.FN%
$include %model_root_path%GAMS/defaults.gms

Positive Variables kapp_slack_ub(j), kapp_slack_lb(j);
Variable kapp_slack_sum;
Equations def_kapp_slack_sum, EnzCap_with_slack(j,j)
; 
def_kapp_slack_sum.. kapp_slack_sum =e= sum(j, kapp_slack_lb(j) + kapp_slack_ub(j));
EnzCap_with_slack(j,j1)$enzload_rxn_coupling(j,j1).. v(j1) * mu =e= (v(j) + kapp_slack_ub(j1) - kapp_slack_lb(j1)) * kapp(j1);

$onMultiR
model rba /all-EnzCap/;
$offMulti
rba.optfile = 1;
solve rba using lp minimizing kapp_slack_sum;

$include %model_root_path%GAMS/default_outputs.gms

file fi_kapp_slack /%gms%.kapp_info.txt/;
put fi_kapp_slack;
fi_kapp_slack.nr=2; fi_kapp_slack.nz=1e-30;
*fi_kapp_slack.lj=2;
put 'rxn':>63 '  ' 'v_ENZLOAD':<14 'kapp_old':<14 'kapp_new':<14 'kapp_slack_ub':<14 'kapp_slack_lb':<14/;
loop(enzload_rxn_coupling(j,j1),
	if ((kapp_slack_ub.l(j1) gt 0) or (kapp_slack_lb.l(j1) gt 0),
		put j1.tl:>63:0 ' ' v.l(j):<14:4 kapp(j1):<14:4 (kapp(j1)*(1+(kapp_slack_ub.l(j1)-kapp_slack_lb.l(j1))/v.l(j))):<14:4 kapp_slack_ub.l(j1):<14:4 kapp_slack_lb.l(j1):<14:4/;
	);
);
putclose fi_kapp_slack;
