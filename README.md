# gambit-publication

Code to reproduce figures and data in {publication name}.


## Instructions

This project is implemented as a [Snakemake](https://snakemake.github.io/) workflow. After
installing and activating the conda environment (see the [Setup](#setup) section below), simply run:

```
snakemake [TARGETS...]
```

from the project's root directory.

By default the `all` target is run, which creates all figures and results. You can also run the
following targets individually:

* `figure_1` through `figure_6`: Generate individual figures.


## Directory structure

* `workflow/`: Snakemake workflow files as well as related Python scripts and Jupyter notebooks.
* `config/`: Workflow configuration files.
* `resources/`: Input data.
  * `genomes/`: Sets of bacterial genomes used for analysis.
  * `gambit-db/`: GAMBIT database files.
* `intermediate-data/`: Output of intermediate workflow targets.
* `results/`: Processed result data.
* `gambit_pub/`: Python package containing common code for this repo.
* `env/`: The conda environment can be installed here.


## Setup

### Required software

All software dependencies are installed using the [conda](https://docs.conda.io) package manager.
If you do not already have it installed, I recommend using the Miniconda installer available
[here](https://docs.conda.io/en/latest/miniconda.html). Make sure the `conda` command is available
in your shell.


### Conda environment

Install the conda environment into the `env/` subdirectory with:

```
conda env create -f env.yaml -p env
```

Before running the workflow you must activate the environment by running `conda activate ./env`
from the project's root directory. This must be done with each new shell session.


### Download source data

Large files in `resources/` are not present in version control and need to be downloaded separately.
You can do this all up front by running the `get_src_data` target. Otherwise the individual data
sets will be downloaded as needed when running the workflow.

