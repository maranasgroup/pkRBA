* `$set solver cplex` recommended if you have CPLEX
$setGlobal solver cplex
$if not setGlobal solver $setGlobal solver soplex
options
	Solver = %solver% /*Solver selection*/
	limrow = 0 /*number of equations listed, 0 is suppresed; turn up if trying to inspect equations (e.g., for infeasibilities) in .lst files*/
	limcol = 0 /*number of variables listed, 0 is suppresed; turn up if trying to inspect variables (e.g., for infeasibilities) in .lst files*/
	iterlim = 1000000 /*iteration limit of solver, for LP it is number of simplex pivots*/
	decimals = 8 /*decimal places for display statement*/
	reslim = 9000 /*wall-clock time limit for solver in seconds*/
	sysout = off /*solver status file report option*/
	solprint = off /*solution printing option*/
* if using cplex, set vmin to the feasibility tolerance of the solver (lowest possible variable value distinguishable from 0), to avoid issues with very small fluxes being treated as infeasible
$if %solver% == cplex $setGlobal vmin 1e-9
* copy solver options file unless keepOptFile is set
$if not setGlobal keepOptFile $call 'rsync -u %system.FP%%solver%.opt .'
