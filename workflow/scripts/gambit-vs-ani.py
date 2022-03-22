"""
Combine pairwise GAMBIT distances and FastANI output into single table.

Expected Snakemake variables:

* input
	gambit: Gambit distance matrix in CSV format.
	fastani: FastANI output.
* output: CSV file to write to.
"""

import ipdb
# ipdb.set_trace()

import numpy as np
import pandas as pd


def ij_to_pair(i, j):
	"""Convert 2D (i, j) index (i > j) to 1D pair index."""
	return i * (i - 1) // 2 + j


### Load data ###

gambit = pd.read_csv(snakemake.input['gambit'], index_col=0)
assert np.array_equal(gambit.index, gambit.columns)

fastani = pd.read_csv(snakemake.input['fastani'], sep='\t', names=['file1', 'file2', 'ani', 'mapped', 'fragments'])


### Make table index ###

ng = gambit.shape[0]
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


### Mean FastANI by pair ###

filename_to_index = {f: i for i, f in enumerate(gambit.index)}
ani1 = np.full(npairs, np.nan)
ani2 = np.full(npairs, np.nan)

for _, row in fastani.iterrows():
	i = filename_to_index[row.file1]
	j = filename_to_index[row.file2]

	if i > j:
		ani1[ij_to_pair(i, j)] = row.ani
	if i < j:
		ani2[ij_to_pair(j, i)] = row.ani

ani_mean = (ani1 + ani2) / 2


### Output table ###

df_out = pd.DataFrame(index=pd.MultiIndex.from_arrays([g1, g2], names=['genome1', 'genome2']))
df_out['gambit'] = gambit.values[g1, g2]
df_out['ani'] = ani_mean

df_out.to_csv(snakemake.output[0])
