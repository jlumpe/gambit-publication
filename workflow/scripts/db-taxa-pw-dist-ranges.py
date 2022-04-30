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

from gambit.db.models import Taxon, only_genomeset
from gambit.db.sqla import file_sessionmaker


### Code ###

def indices_to_slice(indices):
	"""Covert integer arrays of (n ... m) to slice(n, m+1)."""
	if np.array_equal(indices, range(indices[0], indices[-1] + 1)):
		return slice(indices[0], indices[-1] + 1)
	else:
		return indices


### Setup ###

# Get leaf genomes

Session = file_sessionmaker(snakemake.input['db_genomes'])
session = Session()
gset = only_genomeset(session)

leaf_ids = sorted(taxon.id for taxon in gset.taxa.filter(~Taxon.children.any()))
nleaves = len(leaf_ids)
leaf_id_to_index = {id: i for i, id in enumerate(leaf_ids)}

# Genome chunks

chunk_genomes = [pd.read_csv(f) for f in snakemake.input['genome_chunks']]
nchunks = len(chunk_genomes)

chunk_pairs = [(i, j) for i in range(nchunks) for j in range(i, nchunks)]

# Distance matrix chunks

dmat_files = dict()

for file in snakemake.input['dmat_chunks']:
	m = re.search(r'(\d+)-(\d+).h5$', file)
	i = int(m.group(1))
	j = int(m.group(2))
	dmat_files[i, j] = file

assert set(chunk_pairs) == dmat_files.keys()


### Aggregate ###

min_dists = np.full((nleaves, nleaves), np.inf, dtype=np.float32)
max_dists = np.full((nleaves, nleaves), -1, dtype=np.float32)

for c1, c2 in chunk_pairs:
	gb1 = chunk_genomes[c1].groupby('taxon_id')
	gb2 = chunk_genomes[c2].groupby('taxon_id')

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

min_df = pd.DataFrame(min_dists, index=leaf_ids, columns=leaf_ids)
min_df.to_csv(snakemake.output['min_dists'])

max_df = pd.DataFrame(max_dists, index=leaf_ids, columns=leaf_ids)
max_df.to_csv(snakemake.output['max_dists'])
