"""Divide database reference genomes into chunks for distance matrix calculation.

Expected Snakemake variables:

* input
	db_signatures: Database signatures file.
	db_genomes: Database genomes file.
* params
	max_chunk_size: Try to make chunks less than this many genomes.
	taxon_genomes_cap: Cap the number of genomes per leaf taxon to this amount. Used in test mode to
		reduce the amount of work done in computing pairwise distances.
* output
	summary_table: Table of taxon and genome counts by chunk.
	taxa_table: Table of top-level taxa and assignment to chunks.
	chunks_dir: Directory containing genome tables for each chunk.
"""

from pathlib import Path

import pandas as pd

from gambit.db import ReferenceDatabase
from gambit_pub.utils import fix_nullable_int_col


### Setup ###

refdb = ReferenceDatabase.load(snakemake.input['db_genomes'], snakemake.input['db_signatures'])

max_taxon_genomes = snakemake.params['taxon_genomes_cap']

### Table of top-level taxa ###

taxa_rows = []

for taxon in refdb.genomeset.root_taxa():
	taxa_rows.append((
		taxon.id,
		taxon.ncbi_id,
		taxon.name,
		sum(1 for _ in taxon.subtree_genomes()),
	))

taxa_df = pd.DataFrame(taxa_rows, columns=['db_id', 'ncbi_id', 'name', 'ngenomes'])
taxa_df.set_index('db_id', inplace=True)
taxa_df['ncbi_id'] = fix_nullable_int_col(taxa_df['ncbi_id'])


### Assign top-level taxa to chunks ###

taxa_df['chunk'] = None

# Sort by ngenomes then by ID to make sure this is deterministic
taxa_sort = sorted(taxa_df['ngenomes'].items(), key=lambda item: (item[1], item[0]))

current_chunk = 0
current_chunk_size = 0

for db_id, ngenomes in taxa_sort:
	if current_chunk_size + ngenomes > snakemake.params['max_chunk_size']:
		current_chunk += 1
		current_chunk_size = 0

	taxa_df.loc[db_id, 'chunk'] = current_chunk
	current_chunk_size += ngenomes

taxa_df.sort_values(['chunk', 'db_id'])


### Chunks summary table ###

chunk_apply = lambda df: pd.Series(dict(ntaxa=df.shape[0], ngenomes=df['ngenomes'].sum()))
chunks_df = taxa_df.groupby('chunk').apply(chunk_apply)


### Write summary tables ###

chunks_df.to_csv(snakemake.output['summary_table'])
taxa_df.to_csv(snakemake.output['taxa_table'])


### Master genomes table ###

genome_rows = []

for genome, sig_index in zip(refdb.genomes, refdb.sig_indices):
	root_id = genome.taxon.root().id

	genome_rows.append((
		genome.genome_id,
		sig_index,
		taxa_df.loc[root_id, 'chunk'],
		root_id,
		genome.taxon.id,
	))

genomes_df = pd.DataFrame(genome_rows, columns=['id', 'signatures_index', 'chunk', 'root_id', 'taxon_id'])
genomes_df.set_index('id', inplace=True)


### Output per-chunk genome tables ###

chunks_dir = Path(snakemake.output['chunks_dir'])
chunks_dir.mkdir(exist_ok=True)

for chunk, df in genomes_df.groupby('chunk'):

	# Cap number of genomes per taxon
	if max_taxon_genomes is not None:
		df = df.groupby('taxon_id', as_index=False).head(max_taxon_genomes)

	del df['chunk']
	df.to_csv(chunks_dir / f'{chunk}.csv')
