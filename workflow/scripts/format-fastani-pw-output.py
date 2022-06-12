"""
Re-format the output of the FastANI tool to make it a bit easier to work with.

Input is the direct output of the FastANI tool when run with the same set of files as queries and
references (so calculating all pairwise comparisons). The output is a table with one row per unique
pair and an average ANI value (because the result can be slightly different based on which file is
the query and which is the reference). Output index is based on the integer index of the files.

Expected Snakemake variables:

* input: FastANI output.
* params:
  * filenames: Genome file names that were passed to FastANI.
* output: CSV file to write to.
"""

import numpy as np
import pandas as pd


def ij_to_pair(i, j):
	"""Convert 2D (i, j) index (i > j) to 1D pair index."""
	return i * (i - 1) // 2 + j

filenames = snakemake.params['filenames']
ng = len(filenames)
filename_to_index = {f: i for i, f in enumerate(filenames)}


### Load data ###

results = pd.read_csv(snakemake.input[0], sep='\t', names=['file1', 'file2', 'ani', 'mapped', 'fragments'])


### Make table index ###

npairs = ng * (ng - 1) // 2

g1 = np.zeros(npairs, dtype=int)
g2 = np.zeros(npairs, dtype=int)

for i in range(ng):
	for j in range(i):
		p = ij_to_pair(i, j)
		g1[p] = i
		g2[p] = j

# Double-check indexing math is right
assert np.array_equal(list(map(ij_to_pair, g1, g2)), np.arange(npairs))

pair_index = pd.MultiIndex.from_arrays([g1, g2], names=['genome1', 'genome2'])


### Mean FastANI by pair ###

ani1 = np.full(npairs, np.nan)
ani2 = np.full(npairs, np.nan)

for _, row in results.iterrows():
	i = filename_to_index[row.file1]
	j = filename_to_index[row.file2]

	if i > j:
		ani1[ij_to_pair(i, j)] = row.ani
	if i < j:
		ani2[ij_to_pair(j, i)] = row.ani


### Output table ###

out = pd.DataFrame(index=pair_index)
out['ani1'] = ani1
out['ani2'] = ani2
out['ani_mean'] = (ani1 + ani2) / 2

out.to_csv(snakemake.output[0])
