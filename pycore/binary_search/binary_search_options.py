path_gams = '../../GAMS/'
path_pycore = '../../pycore/'
path_enz_mw = '../../GAMS/model/enz_mw_g_per_mmol.txt'
path_pro_mw = '../../GAMS/model/pro_mw_g_per_mmol.txt'

report_file = './binary_search_report.txt' # Text file recording binary search process
mu_tol = 1e-5; # Tolerance of upper and lower bound gap to tolerance search
maxiter = 100; # Maximum number of iteration
mu_min0 = 0; mu_max0 = 0.2; # User-set initial upper and lower bounds of mu
biom_id = 'BIOSYN-BIODILAERO' # Set the ID of the biomass reaction
