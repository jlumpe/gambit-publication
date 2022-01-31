# gambit-publication

Code to reproduce figures and data in {publication name}.


## Overview

This project is based around [GNU Make](https://www.gnu.org/software/make/). Make targets are defined
for all actions such as setting up the environment, downloading source data, and producing specific
figures and results. Rules can be built by running `make NAME_OF_RULE`.
Make keeps track of the dependencies of all targets, so for example building the `results/figure-6`
target will automatically build (TODO - name) to download the required source data if it has not
been done already.

To generate all figures and results in the project you just need to build the `results` target by
running `make results`.


### Make targets

The important targets are:

* `env`: Creates the conda environment containing all software dependencies.
* TODO: download source data
* `results`: Generate all results and figures.

Targets are explained in more detail in the following sections.


### Jupyter notebooks

All analyses are run in Jupyter notebooks in the `notebooks/` directory. All make targets will
execute the required notebooks non-interactively.

To run the notebooks interactively, ... TODO


### Conda environment

All software dependencies, with the exception of those listed in [Prerequisites](#prerequisites),
are installed into a conda environment in the `env/` subdirectory. All make rules use this
environment automatically where necessary, so activating it manually is not required.

To use the environment in your shell for testing/debugging purposes first run `conda activate ./env`
from the project's root directory. This must be done with each new shell session.


### Directory structure

* `notebooks/`: Jupyter notebooks go here.
* `makefiles/`: Sub-makefiles are in this directory, included from the top-level `Makefile`.
* `env/`: The conda environment is installed here.
* `src-data/`: Input data.
  * `genome-sets/`: Sets of bacterial genomes used for analysis.
  * `ncbi/`: Data downloaded from NCBI.
* `results/`: Processed result data.
* `gambit_pub/`: Python package containing common code for this repo.
* `setup/`: Tools to set up the project.


## Setup

### Prerequisites

The following software needs to be installed:

* [GNU Make](https://www.gnu.org/software/make/) - this should already be present in most Linux
  distributions.
* The [conda](https://docs.conda.io) package manager. Recommend to download the Miniconda installer
  [here](https://docs.conda.io/en/latest/miniconda.html).


### Configuration

The `makefiles/env.mk` file contains environment variable definitions that can be customized. You
can edit the file directly or use the `make VAR=VALUE ...` syntax.


### Create the conda environment

The `env` rule creates a conda environment in the `env/` directory with all required software
installed.


### Jupyter notebook kernel

The environment created by this project does not include JupyterLab or the notebook viewer.
Make targets run Jupyter notebooks non-interactively using the `nbconvert` tool. All output
including plots is saved in the notebook file itself, which you can view by opening the file using
a Jupyter installation located in another Conda environment.

In order to run notebooks interactively you will have to add an IPython kernel for this repo's
environment to your Jupyter installation. You can do this by activating the environment and running
the `install-kernel.py` script in the `setup/` directory. Afterwards you can open the notebooks
normally and set the kernel to "GAMBIT Publication".


## Download source data

Large files in `src-data/` are not present in version control and need to be downloaded separately.

TODO make rules


## Generate results

TODO

