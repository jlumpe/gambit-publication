"""
Benchmark performance of GAMBIT vs. FastANI vs. Mash in calculating genomic distance/similarity.

Expected Snakemake variables:
* input
  * queries_fasta_dir
  * queries_list_file
  * refs_fasta_dir
  * refs_list_file
* wildcards
  * genomes_key: Use query/reference genomes defined by config['benchmarks']['dist_benchmarks']['genomes'][genomes_key]
* config: Everything is under the "benchmarks" key.
* output
  * table: CSV file containing main results.
  * extra: JSON file with additional info.
"""

from pathlib import Path
import sys
import json
import shlex

import numpy as np
import pandas as pd
import attrs

import gambit_pub.benchmark as bm
from gambit_pub.utils import truncate_str
from gambit.util.io import read_lines, write_lines


### Code ###

def subsample(n, m, seed=0):
	"""Randomly sample m elements of range(n)."""
	if m is None or m >= n:
		return np.arange(n)

	np.random.seed(seed)
	return np.random.choice(n, m, replace=False)

def log(*args, leading=None, trailing=None, **kw):
	"""Display a log message."""
	if leading is not None:
		print('\n' * leading, file=sys.stderr, end='')
	print(*args, file=sys.stderr, **kw)
	if trailing is not None:
		print('\n' * trailing, file=sys.stderr, end='')

def log_command(cmd, **kw):
	"""Log command before executing it."""
	s = shlex.join(cmd)
	s = truncate_str(s, 400)
	log(f'+ {s}', **kw)


def gambit_signatures_command(list_file, output, params, ncores):
	args = [
		'gambit', 'signatures', 'create',
		'-c', ncores,
		'-k', params['k'],
		'-p', params['prefix'],
		'-l', list_file,
		'-o', output,
		'--no-progress',
	]
	return list(map(str, args))


def gambit_dist_command(query_sigs, ref_sigs, output, ncores):
	args = [
		'gambit', 'dist',
		'-c', ncores,
		'--qs', query_sigs,
		'--rs', ref_sigs,
		'-o', output,
		'--no-progress',
	]
	return list(map(str, args))


def fastani_command(query_list, ref_list, output, params, ncores):
	args = [
		'fastANI',
		'-k', params['k'],
		'--fragLen', params['fraglen'],
		'--threads', ncores,
		'--ql', query_list,
		'--rl', ref_list,
		'-o', output,
	]
	return list(map(str, args))


def mash_sketch_command(files, output, params):
	# Mash sketching does not support parallelization
	args = [
		'mash', 'sketch',
		'-k', params['k'],
		'-s', params['sketch_size'],
		'-o', output,
		*files
	]
	return list(map(str, args))


def mash_dist_command(query_sketch, ref_sketch, params, ncores):
	# Output is recorded in stdout
	args = [
		'mash', 'dist',
		'-p', ncores,
		query_sketch,
		ref_sketch,
	]
	return list(map(str, args))


### Load config ###

conf = snakemake.config['dist_benchmarks']
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
	query_cmd = gambit_signatures_command(query_list_file, qs, params, snakemake.threads)
	log_command(query_cmd)
	dry_run or bm.run(query_cmd)

	rs = ref_sigs[key] = (in_dir / 'refs.h5')
	ref_cmd = gambit_signatures_command(ref_list_file, rs, params, snakemake.threads)
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
	query_cmd = mash_sketch_command(query_files, qs, params)
	log_command(query_cmd)
	dry_run or bm.run(query_cmd)

	rs = ref_sketch[key] = (in_dir / 'refs.msh')
	ref_cmd = mash_sketch_command(ref_files, rs, params)
	log_command(ref_cmd)
	dry_run or bm.run(ref_cmd)


### Commands to be run  ###

runner = bm.BenchmarkRunner(workdir='.')
p = Path('..')

# GAMBIT
for params_key, params in gambit_conf['params'].items():
	for ncores in gambit_conf.get('ncores', default_ncores):
		runner.add_command(
			('gambit', 'query_sigs', params_key, ncores),
			gambit_signatures_command(p / query_list_file, 'query_sigs.h5', params, ncores),
		)
		runner.add_command(
			('gambit', 'ref_sigs', params_key, ncores),
			gambit_signatures_command(p / ref_list_file, 'ref_sigs.h5', params, ncores),
		)
		runner.add_command(
			('gambit', 'dists', params_key, ncores),
			gambit_dist_command(p / query_sigs[params_key], p / ref_sigs[params_key], 'dists.csv', ncores),
		)

# FastANI
for params_key, params in fastani_conf['params'].items():
	for ncores in fastani_conf.get('ncores', default_ncores):
		runner.add_command(
			('fastani', 'fastani', params_key, ncores),
			fastani_command(p / query_list_file, p / ref_list_file, 'out.tsv', params, ncores),
		)

# Mash
for params_key, params in mash_conf['params'].items():
	# Sketching doesn't support multiple cores
	runner.add_command(
		('mash', 'query_sketch', params_key, 1),
		mash_sketch_command([p / f for f in query_files], 'queries.msh', params),
	)
	runner.add_command(
		('mash', 'ref_sketch', params_key, 1),
		mash_sketch_command([p / f for f in ref_files], 'refs.msh', params),
	)

	for ncores in mash_conf.get('ncores', default_ncores):
		runner.add_command(
			('mash', 'dists', params_key, ncores),
			mash_dist_command(p / query_sketch[params_key], p / ref_sketch[params_key], params, ncores),
		)


### Run ####

results = runner.run(replicates=replicates, dry_run=dry_run, shuffle=shuffle)

### Output results ###

if dry_run:
	raise RuntimeError('Intentional rule failure - dry run')

rows = [
	(result.round + 1, *result.key, result.execution_order, *attrs.astuple(result.time))
	for result in results
]
df = pd.DataFrame(rows, columns=['replicate', 'tool', 'command', 'paramset', 'ncores', 'execution_order', 'real', 'user', 'sys'])
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
