"""
To for working with the division of the GAMBIT reference database genomes into multiple chunks, as
used in the the `db-pw-dists.smk` rules file.
"""

import os
import re
from contextlib import contextmanager

import numpy as np
import pandas as pd
import h5py as h5

from gambit.db import ReferenceGenomeSet, Taxon, load_genomeset


def _get_taxon(taxon_or_tid):
	return taxon_or_tid if isinstance(taxon_or_tid, Taxon) else taxa_by_id[taxon_or_tid]

def _get_taxon_id(taxon_or_tid):
	return taxon_or_tid.id if isinstance(taxon_or_tid, Taxon) else taxon_or_tid


class GenomeChunks:
	"""Division of database genomes into multiple chunks.

	This basically represents the full output of the db_pw_dists_make_chunks checkpoint/rule.

	Attributes
	----------
	n : int
		Number of chunks
	pairs
		All ordered pairs of chunk indices.
	taxa_df : pandas.DataFrame
		taxa_table output of rule.
	genome_dfs : List[pandas.DataFrame]
		Tables in chunks_dir directory output of rule (one table per chunk).
	genomeset : gambit.db.models.ReferenceGenomeSet
		Optional, enables lookups based on taxonomy tree structure.
	"""

	def __init__(self, taxa_df, genome_dfs, genomeset=None):
		self.taxa_df = taxa_df
		self.genome_dfs = list(genome_dfs)
		self.n = len(self.genome_dfs)
		self.pairs = tuple((i, j) for i in range(self.n) for j in range(i, self.n))

		self.genomeset = genomeset
		if genomeset is not None:
			self.taxa = genomeset.taxa
			self.root_taxa = genomeset.root_taxa()
			self.taxa_by_id = {taxon.id: taxon for taxon in genomeset.taxa}
		else:
			self.taxa = None
			self.root_taxa = None
			self.taxa_by_id = None

	@classmethod
	def load(cls, taxa_file, genome_dir, genomes_db_file=None):
		"""Load from csv files."""

		taxa_df = pd.read_csv(taxa_file, index_col=0)
		n = taxa_df['chunk'].max() + 1

		genome_files = [os.path.join(genome_dir, f'{i}.csv') for i in range(n)]
		genome_dfs = [pd.read_csv(f) for f in genome_files]

		if genomes_db_file is None:
			gset = None
		else:
			_session, gset = load_genomeset(genomes_db_file)

		return cls(taxa_df, genome_dfs, gset)

	def get_taxon_chunk(taxon):
		return chunk_taxa.loc[_get_taxon(taxon).root().id, 'chunk']

	def get_genome_indices_direct(chunk, taxa_or_tids, asint=True):
		bools = np.in1d(chunk_genomes.loc[chunk, 'taxon_id'], tids)
		return np.flatnonzero(bools) if asint else bools

	def get_genome_indices_subtree(chunk, taxon, asint=True):
		subtree_tids = [t.id for t in taxon.traverse()]
		return get_genome_indices_direct(chunk, subtree_tids, asint)


def map_dmat_chunk_files(paths):
	"""Decode the chunk pair for output files of the db_pw_dists_chunk rule based on file names."""
	d = dict()

	for path in paths:
		fname = os.path.basename(path).split('.')[0]
		m = re.match(r'(\d+)-(\d+)', fname)
		if m is None:
			raise ValueError(f'File name does not have the expected format: {fname}')
		c1, c2 = map(int, m.groups())
		assert c1 <= c2
		d[c1, c2] = path

	return d


class ChunkedDists:
	"""
	Attributes
	----------
	chunks : GenomeChunks
	dmat_files : Dict[Tuple, str]
		Mapping from chunk index pairs to files containing distance matrix chunks.
	"""

	def __init__(self, chunks: GenomeChunks, dmat_files):
		self.chunks = chunks
		self.dmat_files = dmat_files

	# def open_dmat(self, i, j):
