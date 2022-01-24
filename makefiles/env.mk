# Customizable environment variables
# Override by editing this file or from the command line like "make VAR=VALUE ..."

# Set this to a non-empty value to enable some options that help in developing this repo
DEV=

# Conda command or path to conda executable
CONDA=conda


# Don't edit variables below:

# Run a command in the conda environment
CONDA_RUN=$(CONDA) run -p ./env
# Executes a Jupyter notebook noninteractively
RUN_NOTEBOOK=jupyter nbconvert --to notebook --execute --inplace
