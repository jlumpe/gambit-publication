# GAMBIT (Genomic Approximation Method for Bacterial Identification and Tracking): A methodology to rapidly leverage whole genome sequencing of bacterial isolates for clinical identification

A [Snakemake](https://snakemake.github.io/) workflow to (re)produce figures and data in the initial
GAMBIT [publication](https://doi.org/10.1101/2022.06.14.496173). Source code for GAMBIT itself is
located [here](https://github.com/hesslab-gambit/gambit/).


## Instructions

After  installing and activating the conda environment (see the [Setup](#setup) section below),
simply run:

```
snakemake [TARGETS...]
```

from the project's root directory. `TARGETS` are one or more rule names or output files. By default
the `main` rule is run, which creates figures 1-6. See the [Targets](#targets) section for a list
of options.


## Directory structure

* `workflow/`: Snakemake workflow files and related scripts.
* `config/`: Workflow configuration files.
* `resources/`: Input data.
  * `genomes/`: Sets of bacterial genomes used for analysis.
  * `gambit-db/`: GAMBIT database files.
* `intermediate-data/`: Output of intermediate workflow targets.
* `results/`: Processed result data.
* `gambit_pub/`: Python package containing common code for this repo.
* `env/`: The conda environment can be installed here.


## Setup

This workflow has been built and tested for Linux only. It may work on Mac (haven't tested) but I
believe there are issues preventing it from running on Windows.

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

### Install GAMBIT

The preferred way to install GAMBIT is through the Bioconda channel:

```
conda install -c bioconda hesslab-gambit
```

Make sure your Conda environment is activated first.

### Configuration ###

Most editable config settings are in `config/config.yaml`.

### Download source data

Large files in `resources/` are not present in version control and need to be downloaded separately.
You can do this all up front by running the `fetch_src_data` target, which may make things easier to
debug if you run into any connection problems. Otherwise the individual data sets will be
downloaded as needed when running the workflow.


## Targets

This is a list of all "endpoint" rules and output files which you may want to run.
It does not include rules which generate intermediate data.

### Aggregate rules

| Rule             | Description                                                                    |
|------------------|--------------------------------------------------------------------------------|
| `all`            | `main` and `supplemental`.                                                     |
| `main`           | Generate all primary figures.                                                  |
| `supplemental`   | Generate all supplemental figures. Note - suppplemental figure 1 is VERY slow. |
| `fetch_src_data` | Download all source data.                                                      |

### Main results

| Rule       | Output                               | Description         |
|------------|--------------------------------------|---------------------|
| `figure_1` | `results/figures/figure-1.{png,csv}` | Generate figure 1.  |
| `figure_2` | `results/figures/figure-2{a,b}.png`  | Generate figure 2.  |
| `figure_3` | `results/figures/figure-3.png`       | Generate figure 3.  |
| `figure_4` | `results/figures/figure-4{a,b}.png`  | Generate figure 4.  |
| `figure_5` | `results/figures/figure-5{a,b}.png`  | Generate figure 5.  |
| `figure_6` | `results/figures/figure-6.png`       | Generate figure 6.  |

### Supplemental results

| Rule                     | Output                                      | Description                                       |
|--------------------------|---------------------------------------------|---------------------------------------------------|
| `supplemental_figure_1`  | `results/figures/supplemental-figure-1.png` | Generate supplemental figure 1. Note - VERY slow. |
| `supplemental_figure_2`  | `results/fibures/supplemental-figure-2.png` | Generate supplemental figure 2.                   |

### Benchmarks

| Rule              | Output                             | Description                                         |
|-------------------|------------------------------------|-----------------------------------------------------|
| `benchmark_query` | `results/benchmarks/gambit-query/` | Benchmark GAMBIT taxonomic classification from CLI. |

### Source data

| Rule                 | Output                              | Description                                                                     |
|----------------------|-------------------------------------|---------------------------------------------------------------------------------|
| `fetch_gambit_db`    | `resources/gambit-db/`              | Download GAMBIT reference database files.                                       |
|                      | `resources/genomes/set{1,2}/fasta/` | Download FASTA files for data set 1 or 2 from NCBI. Invoke by output directory. |
|                      | `resources/genomes/set{3,4}/fasta/` | Download FASTA files for data set 3 or 4. Invoke by output directory.           |
| `fetch_genome_set_5` | `resources/genomes/set5/fasta/`     | Download FASTA files for data set 5.                                            |


## Development

You can enable "test mode" by adding `--config test=1` to the command line options. This loads
an alternate set of parameters which greatly reduces the amount of work to be done.
