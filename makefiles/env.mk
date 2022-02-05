# Customizable environment variables
# Override by editing this file or from the command line like "make VAR=VALUE ..."

# Set this to a non-empty value to enable some options that help in developing this repo
devmode=

# Conda command or path to conda executable
conda=conda


# Don't edit variables below:

# Run a command in the conda environment
conda_run = $(conda) run -p $(CURDIR)/env
# Executes a Jupyter notebook non-interactively. Need to prepend with $(conda_run).
run_notebook = jupyter nbconvert --to notebook --execute --inplace


# Download a private file from GCS. The environment variable GCS_OAUTH2_TOKEN must be set.
# This is temporary until we get everything set up on the public GCS bucket
# Use as $(call DOWNLOAD_GCS,object[,destination])
DOWNLOAD_GCS = wget --header "Authorization: Bearer $(GCS_OAUTH2_TOKEN)" \
	"https://storage.googleapis.com/$(1)" $(if $2,-O "$2")
