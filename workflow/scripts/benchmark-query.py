"""
Benchmark performance of GAMBIT query command.

Expected Snakemake variables:
* input
  * db_genomes: GAMBIT database genomes file.
  * db_signatures: GAMBIT database signatures file.
  * query_file_dirs: Directory containing query files for each genome set.
  * query_list_files: File containing list of query file names for each genome set.
* params
  * genome_sets: Labels for genome sets used.
* output
  * table: CSV file containing main results.
  * extra: JSON file with additional info.
"""

from pathlib import Path

import pandas as pd
import attrs

import gambit_pub.benchmark as bm
from gambit.util.misc import zip_strict


### Setup ###

db_dir = Path(snakemake.input['db_genomes']).parent.absolute()
fasta_dirs = snakemake.input['query_file_dirs']
list_files = snakemake.input['query_list_files']

conf = snakemake.config['query_benchmark']
ncores = conf['ncores']
replicates = conf['replicates']
shuffle = conf['shuffle']

dry_run = bool(snakemake.config['benchmark_dry_run'])

genome_sets = snakemake.params['genome_sets']


### Run ###

runner = bm.BenchmarkRunner(workdir='.')
p = Path('..')

for nc in ncores:
	for gset, fasta_dir, list_file in zip_strict(genome_sets, fasta_dirs, list_files):
		args = [
			'gambit',
			'--db', db_dir,
			'query',
			'-l', p / list_file,
			'--ldir', p / fasta_dir,
			'--no-progress',
		]
		if nc is not None:
			args.extend(['-c', int(nc)])
		runner.add_command((gset, nc), args)

results = runner.run(replicates=replicates, dry_run=dry_run, shuffle=shuffle)


### Output results ###

if dry_run:
	raise RuntimeError('Intentional rule failure - dry run')

rows = [
	(result.round + 1, *result.key, result.execution_order, *attrs.astuple(result.time))
	for result in results
]
df = pd.DataFrame(rows, columns=['replicate', 'genome_set', 'ncores', 'execution_order', 'real', 'user', 'sys'])
df.sort_values(['replicate', 'genome_set', 'ncores'], inplace=True)
df.to_csv(snakemake.output[0], index=False)
