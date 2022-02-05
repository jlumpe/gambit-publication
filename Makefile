# File name of this makefile
# Other makefiles can check if this is defined to see if they were included properly
ROOT_MAKEFILE := $(lastword $(MAKEFILE_LIST))

SHELL = /bin/bash


include makefiles/env.mk
include makefiles/src-data.mk
include makefiles/results.mk


.PHONY: list env


# List available targets
# From https://stackoverflow.com/a/26339924/1775059
list:
	@LC_ALL=C $(MAKE) -pRrq -f $(ROOT_MAKEFILE) : 2>/dev/null | awk -v RS= -F: '/^# File/,/^# Finished Make data base/ {if ($$1 !~ "^[#.]") {print $$1}}' | sort | egrep -v -e '^[^[:alnum:]]' -e '^$@$$'


# Create the conda environment
env:
	$(conda) env create -f env.yaml --prefix ./env
	$(conda_run) pip $(if $(devmode), "-e") install ./gambit_pub
