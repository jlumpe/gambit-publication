"""Calculate minimum and maximum genome distances aggregated by leaf taxon.

Expected Snakemake variables:

* input
	db_genomes: Database genomes file.
	chunks_dir: Directory containing genome tables for each chunk.
* output
"""

import re

import numpy as np
import pandas as pd
import h5py as h5

from gambit.db import Taxon, load_genomeset
from gambit_pub.db_pw_dists import GenomeChunks, map_dmat_chunk_files


### Code ###

def indices_to_slice(indices):
	"""Covert integer arrays of (n ... m) to slice(n, m+1)."""
	if np.array_equal(indices, range(indices[0], indices[-1] + 1)):
		return slice(indices[0], indices[-1] + 1)
	else:
		return indices


### Setup ###

# Get leaf genomes

session, gset = load_genomeset(snakemake.input['db_genomes'])

leaf_ids = sorted(taxon.id for taxon in gset.taxa.filter(~Taxon.children.any()))
nleaves = len(leaf_ids)
leaf_id_to_index = {id: i for i, id in enumerate(leaf_ids)}


# Genome chunks

chunks = GenomeChunks.load(snakemake.input['chunk_taxa'], snakemake.input['chunk_genomes_dir'])

dmat_files = map_dmat_chunk_files(snakemake.input['dmat_chunks'])
assert dmat_files.keys() == set(chunks.pairs)


### Aggregate ###

min_dists = np.full((nleaves, nleaves), np.inf, dtype=np.float32)
max_dists = np.full((nleaves, nleaves), -1, dtype=np.float32)

from tqdm import tqdm
for c1, c2 in tqdm(chunks.pairs):
	gb1 = chunks.genome_dfs[c1].groupby('taxon_id')
	gb2 = chunks.genome_dfs[c2].groupby('taxon_id')

	with h5.File(dmat_files[c1, c2]) as file:
		dmat = file['dmat'][:]

	for tid1, indices1 in gb1.indices.items():
		li1 = leaf_id_to_index[tid1]
		indices1 = indices_to_slice(indices1)

		for tid2, indices2 in gb2.indices.items():
			li2 = leaf_id_to_index[tid2]
			indices2 = indices_to_slice(indices2)

			subdmat = dmat[indices1, indices2]
			min_dists[li1, li2] = min_dists[li2, li1] = min(min_dists[li1, li2], subdmat.min())
			max_dists[li1, li2] = max_dists[li2, li1] = max(max_dists[li1, li2], subdmat.max())


### Write output ###

minmax_index = pd.Series(leaf_ids, name='taxon_id')

min_df = pd.DataFrame(min_dists, index=minmax_index, columns=leaf_ids)
min_df.to_csv(snakemake.output['min_dists'])

max_df = pd.DataFrame(max_dists, index=minmax_index, columns=leaf_ids)
max_df.to_csv(snakemake.output['max_dists'])
