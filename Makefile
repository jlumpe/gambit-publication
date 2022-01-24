include makefiles/env.mk
include makefiles/src-data.mk
include makefiles/results.mk


SHELL=/bin/bash


.PHONY: env


# Create the conda environment
env:
	$(CONDA) env create -f env.yaml --prefix ./env
	$(CONDA_RUN) pip $(if $(DEV), "-e") install ./gambit_pub
