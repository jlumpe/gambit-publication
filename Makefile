include makefiles/env.mk
include makefiles/src-data.mk
include makefiles/results.mk


SHELL=/bin/bash

# Run a command in the conda environment
CONDA_RUN=$(CONDA) run -p ./env
# Executes a Jupyter notebook noninteractively
RUN_NOTEBOOK=jupyter nbconvert --to notebook --execute --inplace


.PHONY: env


# Create the conda environment
env:
	$(CONDA) env create -f env.yaml --prefix ./env
	$(CONDA_RUN) pip $(if $(DEV), "-e") install ./gambit_pub
