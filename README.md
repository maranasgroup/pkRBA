# pkRBA
Resource Balance Analysis model for *Pichia kudriavzevii* SD108

The updated genome-scale metabolic model used to build pkRBA is available [here](build_model/input/non-decomp/iPickudri866.json).
# Installation
To install directly via the environment.yml file, Anaconda must already be installed. To create a self-contained environment for RBA, activate the Anaconda prompt and type:
```conda env create -f <path to environment.yml file>```
# Picking a solver
By default, CPLEX is used as the solver for running the RBA model. To utilize a different solver (e.g., if you do not have a CPLEX license), edit [GAMS/default_solver.gms](GAMS/default_solver.gms).
# Building an RBA model
To create a new RBA model or update the one from this repository, see [build_model/README.md](build_model/README.md).
# Running yield calculations
- After building a new RBA model, update the following files: 
  - [application/input/prod_mw.xlsx](application/input/prod_mw.xlsx): fill out `product` and `MW (g/mol)` columns. All others are optional.
  - [application/input/pathways.xlsx](application/input/pathways.xlsx): if any reactions must be added to allow making your product that are heterologous (i.e., not normally present) add them here (including which products they can be used for, the full rxn ID, and the stoichiometry). Best to avoid including any reactions that are already present in your model; can run [application/A1_run_yield_predictions.py](application/A1_run_yield_predictions.py) to check for potentially redundant ones.
  - [application/input/phenotype_common.txt](application/input/phenotype_common.txt): If you have any conditions you'd like to enforce during yield calculations for all products, enforce them here.
- To run yields:
  - Go under [application/](application/) and the folder for the type of yield calculations you want to run (e.g., `output_max_default` for max yields using default RBA settings). Run your chosen folder's `run.sh` file.
# Funding
This work was funded by the DOE Center for Advanced Bioenergy and Bioproducts Innovation (U.S. Department of Energy, Office of Science, Biological and Environmental Research Program under Award Number DE-SC0018420). Any opinions, findings, and conclusions or recommendations expressed in this publication are those of the author(s) and do not necessarily reflect the views of the U.S. Department of Energy. Computations for this research were performed on the Pennsylvania State University’s Roar Collab supercomputer. The authors of this work recognize the Penn State Institute for Computational and Data Sciences (RRID:SCR_025154) for providing access to computational research infrastructure within the Roar Core Facility (RRID: SCR_026424).
