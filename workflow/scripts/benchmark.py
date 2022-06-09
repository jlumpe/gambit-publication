"""
Benchmark performance of GAMBIT vs. FastANI vs. Mash.

Expected Snakemake variables:
* input
  * queries_fasta_dir
  * queries_list_file
  * refs_fasta_dir
  * refs_list_file
* wildcards
  * genomes_key: Use query/reference genomes defined by config['benchmarks']['genomes'][genomes_key]
* config: Everything is under the "benchmarks" key.
* output
  * table: CSV file containing main results.
  * extra: JSON file with additional info.
"""

from pathlib import Path
import sys
from subprocess import CalledProcessError
import json

import numpy as np
import pandas as pd
import attrs

import gambit_pub.benchmark as bm
from gambit.util.io import read_lines, write_lines


### Code ###

def subsample(n, m, seed=0):
	"""Randomly sample m elements of range(n)."""
	if m is None or m >= n:
		return np.arange(n)

	np.random.seed(seed)
	return np.random.choice(n, m, replace=False)

def truncate_str(s, n):
	l = len(s)
	if l <= n:
		return s
	return f'{s[:n]}... ({l - n} characters omitted)'

def log(*args, leading=None, trailing=None, **kw):
	"""Display a log message."""
	if leading is not None:
		print('\n' * leading, file=sys.stderr, end='')
	print(*args, file=sys.stderr, **kw)
	if trailing is not None:
		print('\n' * trailing, file=sys.stderr, end='')

def log_command(cmd: bm.Command, **kw):
	"""Log command before executing it."""
	s = truncate_str(str(cmd), 400)
	log(f'+ {s}', **kw)


### Load config ###

conf = snakemake.config['benchmarks']
default_ncores = conf['ncores']
replicates = conf['replicates']
shuffle = conf['shuffle']

genome_params = conf['genomes'][snakemake.wildcards.genomes_key]
queries_subsample = genome_params['queries']['subsample']
refs_subsample = genome_params['refs']['subsample']

gambit_conf = conf['gambit']
fastani_conf = conf['fastani']
mash_conf = conf['mash']

dry_run = bool(snakemake.config['benchmark_dry_run'])


### Locate genome files ###

queries_dir = Path(snakemake.input.queries_fasta_dir)
query_filenames = np.asarray(list(read_lines(snakemake.input.queries_list_file)))
query_indices = subsample(len(query_filenames), queries_subsample)
query_filenames = query_filenames[query_indices]
query_files = [queries_dir / f for f in query_filenames]

refs_dir = Path(snakemake.input.refs_fasta_dir)
ref_filenames = np.asarray(list(read_lines(snakemake.input.refs_list_file)))
ref_indices = subsample(len(ref_filenames), refs_subsample)
ref_filenames = ref_filenames[ref_indices]
ref_files = [refs_dir / f for f in ref_filenames]


### Create input files for commands ###
# Generate the inputs once at the start so we can run the commands in any order

log('Running benchmarks in temporary directory', Path('.').absolute())

# List files
query_list_file = Path('queries.txt')
write_lines(query_files, query_list_file)

ref_list_file = Path('refs.txt')
write_lines(ref_files, ref_list_file)

# GAMBIT signatures
query_sigs = dict()
ref_sigs = dict()
for key, params in gambit_conf['params'].items():
	log(f'Generating input files for GAMBIT parameter group {key!r}:', leading=1)
	in_dir = Path(f'gambit-input-{key}')
	dry_run or in_dir.mkdir()

	qs = query_sigs[key] = (in_dir / 'queries.h5')
	query_cmd = bm.gambit_signatures_command(query_list_file, qs, params, snakemake.threads)
	log_command(query_cmd)
	dry_run or bm.run(query_cmd)

	rs = ref_sigs[key] = (in_dir / 'refs.h5')
	ref_cmd = bm.gambit_signatures_command(ref_list_file, rs, params, snakemake.threads)
	log_command(ref_cmd)
	dry_run or bm.run(ref_cmd)

# Mash sketches
query_sketch = dict()
ref_sketch = dict()
for key, params in mash_conf['params'].items():
	log(f'Generating input files for Mash parameter group {key!r}:', leading=1)
	in_dir = Path(f'mash-input-{key}')
	dry_run or in_dir.mkdir()

	qs = query_sketch[key] = (in_dir / 'queries.msh')
	query_cmd = bm.mash_sketch_command(query_files, qs, params)
	log_command(query_cmd)
	dry_run or bm.run(query_cmd)

	rs = ref_sketch[key] = (in_dir / 'refs.msh')
	ref_cmd = bm.mash_sketch_command(ref_files, rs, params)
	log_command(ref_cmd)
	dry_run or bm.run(ref_cmd)


### Commands to be run  ###

commands = dict()
p = Path('..')

# GAMBIT
for params_key, params in gambit_conf['params'].items():
	for ncores in gambit_conf.get('ncores', default_ncores):
		commands['gambit', 'query_sigs', params_key, ncores] = \
			bm.gambit_signatures_command(p / query_list_file, 'query_sigs.h5', params, ncores)
		commands['gambit', 'ref_sigs', params_key, ncores] = \
			bm.gambit_signatures_command(p / ref_list_file, 'ref_sigs.h5', params, ncores)
		commands['gambit', 'dists', params_key, ncores] = \
			bm.gambit_dist_command(p / query_sigs[params_key], p / ref_sigs[params_key], 'dists.csv', ncores)

# FastANI
for params_key, params in fastani_conf['params'].items():
	for ncores in fastani_conf.get('ncores', default_ncores):
		commands['fastani', 'fastani', params_key, ncores] = \
			bm.fastani_command(p / query_list_file, p / ref_list_file, 'out.tsv', params, ncores)

# Mash
for params_key, params in mash_conf['params'].items():
	# Sketching doesn't support multiple cores
	commands['mash', 'query_sketch', params_key, 1] = \
		bm.mash_sketch_command([p / f for f in query_files], 'queries.msh', params)
	commands['mash', 'ref_sketch', params_key, 1] = \
		bm.mash_sketch_command([p / f for f in ref_files], 'refs.msh', params)

	for ncores in mash_conf.get('ncores', default_ncores):
		commands['mash', 'dists', params_key, ncores] = \
			bm.mash_dist_command(p / query_sketch[params_key], p / ref_sketch[params_key], params, ncores)


### Run ####

items = []

for r in range(replicates):
	log(f'*** Starting round {r + 1} of {replicates} ***', leading=2)

	# Shuffle list of commands
	commands_iter = np.asarray(list(commands.items()), dtype=object)
	if shuffle:
		np.random.seed(r)
		np.random.shuffle(commands_iter)

	for i, (key, command) in enumerate(commands_iter):
		keystr = '-'.join(map(str, key))
		log(f'Running {keystr} ({i + 1} of {len(commands)})', leading=1)
		log_command(command)

		if dry_run:
			continue

		workdir = Path(f'{keystr}-{r}')
		workdir.mkdir()

		# Try running
		try:
			result = bm.benchmark_command(
				command, workdir, 'time',
				stdout=workdir / 'stdout',
				stderr=workdir / 'stderr',
			)

		except CalledProcessError as e:
			# Print stdout and stderr and re-raise
			with open(workdir / 'stdout') as f:
				log(f.read())
			with open(workdir / 'stderr') as f:
				log(f.read())
			raise

		log(f'Completed in {result.time}')

		items.append((r, *key, i, result))


### Output results ###

if dry_run:
	raise RuntimeError('Intentional rule failure - dry run')

rows = [
	(r + 1, *key, *attrs.astuple(result.time))
	for r, *key, result in items
]
df = pd.DataFrame(rows, columns=['replicate', 'tool', 'command', 'paramset', 'ncores', 'run_order', 'real', 'user', 'sys'])
df.sort_values(['replicate', 'tool', 'command', 'paramset', 'ncores'], inplace=True)
df.to_csv(snakemake.output['table'], index=False)

extra = dict(
	query_indices=list(map(int, query_indices)),
	query_files=list(query_filenames),
	ref_indices=list(map(int, ref_indices)),
	ref_files=list(ref_filenames),
)
with open(snakemake.output['extra'], 'w') as f:
	json.dump(extra, f)
