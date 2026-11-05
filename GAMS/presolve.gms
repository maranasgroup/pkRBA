* common routines for preparing for solving (e.g., ways of removing blocked rxns)

$include %model_root_path%GAMS/find_blocked_rxns.gms

* for rxns where v_sm_lo() is set to 0, prevent their reverse directions from being used (to improve solving time and avoid infeasibilities), and same w/ v_sm_up(); turns off if simulating cycling, but may be worth turning off in other cases too (e.g., this line sometimes causes bugs when v_sm values are updated later on)
v.fx(j)$(sum(sm_j$(v_sm_lo(sm_j) ge 0 and dir(sm_j, j) = -1 and cap_gross_flux(sm_j) ne 0), 1)) = 0;
v.fx(j)$(sum(sm_j$(v_sm_up(sm_j) le 0 and dir(sm_j, j) = 1 and cap_gross_flux(sm_j) ne 0), 1)) = 0;

* prevent enzloads from being used if their rxns are blocked from carrying flux, unless they have to be used (i.e., if v.lo != 0)
v.fx(j)$(enzload(j) and v.lo(j) eq 0 and smax(j1$enzload_rxn_coupling(j1, j), v.up(j1)) eq 0) = 0;
