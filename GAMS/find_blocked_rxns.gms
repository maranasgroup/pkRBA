* finds several kinds of blocked rxns; must be loaded in by another .gms file

scalar nonzero_ub, rxns_blocked_count /1/;
set rxns_with_met(j);

while(rxns_blocked_count gt 0,
	rxns_blocked_count = 0;
* find metabolites that can't be used, and block all their associated rxns
	loop(i,
		rxns_with_met(j)=no;
		rxns_with_met(j)$(S(i,j) and v.up(j) gt %vmin%)=yes;
*	put log; put i.tl:0/;
*	loop(j$rxns_with_met(j),
*		put system.tab, j.tl:0/;
*	);
*	putclose;
		if(card(rxns_with_met) eq 0,
			continue;
		);
		if(card(rxns_with_met) le 1,
			rxns_blocked_count = rxns_blocked_count + 1;
			v.fx(rxns_with_met) = 0;
*		put log; put i.tl:0/; putclose;
		elseif(smax(j$rxns_with_met(j),S(i,j)) le 0 or smin(j$rxns_with_met(j),S(i,j)) ge 0),
			rxns_blocked_count = rxns_blocked_count + 1;
			v.fx(rxns_with_met) = 0;
		);
	);


* if all rxns in sm_j are blocked, block it too (and vice versa)
	loop(sm_j,
		if((v_sm_up(sm_j) le %vmin%) and (v_sm_lo(sm_j) ge -%vmin%),
			loop(j$(dir(sm_j,j) ne 0 and v.up(j) gt %vmin%),
				v.fx(j) = 0;
				rxns_blocked_count = rxns_blocked_count + 1;
			);
		else
			if(cap_gross_flux(sm_j) ne 0,
				if(v_sm_lo(sm_j) ge -%vmin%,
					rxns_blocked_count = rxns_blocked_count + sum(j$(dir(sm_j,j) = -1 and v.up(j) gt %vmin%),1);
					v.fx(j)$(dir(sm_j,j) = -1) = 0;
				elseif(v_sm_up(sm_j) le %vmin%),
					rxns_blocked_count = rxns_blocked_count + sum(j$(dir(sm_j,j) = 1 and v.up(j) gt %vmin%),1);
					v.fx(j)$(dir(sm_j,j) = 1) = 0;
				);
			);
			nonzero_ub = 0;
			if(smax(j$(dir(sm_j,j) ne 0),v.up(j)) gt %vmin%,
				nonzero_ub = 1;
			);
			if(nonzero_ub = 0,
				v_sm_lo(sm_j) = 0;
				v_sm_up(sm_j) = 0;
				rxns_blocked_count = rxns_blocked_count + 1;
				put log; put 'blocked ' sm_j.tl:0 ' due to lack of v(j) entries with bounds over ' %vmin%; putclose;
			);
		);
	);
);
