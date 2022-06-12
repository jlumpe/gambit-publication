"""
Combine pairwise GAMBIT distances and FastANI output into single table.

Expected Snakemake variables:

* input
  * gambit: Gambit distance matrix in CSV format.
  * fastani: Formatted FastANI output from format_fastani_results rule.
* output:
  * pairs: CSV file of GAMBIT and FastANI scores for all genome pairs.
  * stats: JSON file of statistics, including Pearson and Spearman correlation.
"""

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

import gambit.util.json as gjson


### Load data ###

gambit = pd.read_csv(snakemake.input['gambit'], index_col=0)
assert np.array_equal(gambit.index, gambit.columns)

fastani = pd.read_csv(snakemake.input['fastani'], index_col=[0, 1])

ng = gambit.shape[0]
npairs = ng * (ng - 1) // 2
assert fastani.shape[0] == npairs


### Combine ###

combined = fastani['ani_mean'].to_frame('ani')
combined['gambit'] = np.fromiter(
	(gambit.values[i, j] for i, j in combined.index),
	dtype=np.float32,
	count=npairs
)


### Statistics ###

both_reported = ~pd.isnull(fastani['ani1']) & ~pd.isnull(fastani['ani2'])
either_reported = ~pd.isnull(fastani['ani1']) | ~pd.isnull(fastani['ani2'])
x = combined['gambit'][both_reported]
y = combined['ani'][both_reported]

stats = dict(
	ngenomes=ng,
	npairs=npairs,
	ani_both_reported=both_reported.sum(),
	ani_either_reported=either_reported.sum(),
	pearson=pearsonr(x, y)[0],
	spearman=spearmanr(x, y)[0],
)


### Output ###

combined.to_csv(snakemake.output['pairs'])

with open(snakemake.output['stats'], 'wt') as f:
	gjson.dump(stats, f)
