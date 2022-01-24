# Customizable environment variables
# Override by editing this file or from the command line like "make VAR=VALUE ..."

# Set this to a non-empty value to enable some options that help in developing this repo
DEV=

# Conda command or path to conda executable
CONDA=conda

# URLs to download GAMBIT database
GAMBIT_DB_BASE_URL=https://storage.googleapis.com/hesslab-gambit-public/databases/refseq-curated/1.0-beta
GAMBIT_DB_GENOMES_URL=$(GAMBIT_DB_BASE_URL)/gambit-genomes-1.0b1-210719.db
GAMBIT_DB_SIGNATURES_URL=$(GAMBIT_DB_BASE_URL)/gambit-signatures-1.0b1-210719.h5


# Don't edit variables below:

# Run a command in the conda environment
CONDA_RUN=$(CONDA) run -p ./env
# Executes a Jupyter notebook noninteractively
RUN_NOTEBOOK=jupyter nbconvert --to notebook --execute --inplace
