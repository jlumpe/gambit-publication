"""Divide database reference genomes into chunks for distance matrix calculation.

Expected Snakemake variables:

* input
	db_signatures: Database signatures file.
	db_genomes: Database genomes file.
* params
	max_chunk_size: Try to make chunks less than this many genomes.
* output
	summary_table: Table of taxon and genome counts by chunk.
	taxa_table: Table of top-level taxa and assignment to chunks.
	chunks_dir: Directory containing genome tables for each chunk.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from gambit.db import ReferenceDatabase


def subtree_genomes(taxon):
	return [genome for t in taxon.traverse() for genome in t.genomes]

def ancestor_of_rank(taxon, rank: str):
	"""Get ID of taxon's ancestor with given rank."""
	for ancestor in taxon.ancestors(True):
		if ancestor.rank == rank:
			return ancestor.id
	return None

def fix_nullable_int_col(values):
	"""Fix column containing ints/Nones after Pandas coerces it to float data type."""
	return np.asarray([None if pd.isnull(v) else int(v) for v in values], dtype=object)


### Open database ###

refdb = ReferenceDatabase.load(snakemake.input['db_genomes'], snakemake.input['db_signatures'])


### Table of top-level taxa ###

taxa_rows = []

for taxon in refdb.genomeset.root_taxa():
	taxa_rows.append((
		taxon.id,
		taxon.ncbi_id,
		taxon.name,
		len(subtree_genomes(taxon))
	))

taxa_df = pd.DataFrame(taxa_rows, columns=['db_id', 'ncbi_id', 'name', 'ngenomes'])
taxa_df.set_index('db_id', inplace=True)
taxa_df['ncbi_id'] = fix_nullable_int_col(taxa_df['ncbi_id'])
taxa_df.sort_values('ngenomes', inplace=True)


### Assign top-level taxa to chunks ###

taxa_df['chunk'] = None

current_chunk = 0
current_chunk_size = 0

for db_id, ngenomes in taxa_df['ngenomes'].items():
	if current_chunk_size + ngenomes > snakemake.params['max_chunk_size']:
		current_chunk += 1
		current_chunk_size = 0

	taxa_df.loc[db_id, 'chunk'] = current_chunk
	current_chunk_size += ngenomes


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
	del df['chunk']
	df.to_csv(chunks_dir / f'{chunk}.csv')
