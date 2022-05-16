"""
Generate table of GAMBIT distance/ANI correlations.

Expected Snakemake variables:

* input: gambit_vs_ani rule output for each row of parameters table.
* params
  * params_df: Data frame of parameter combinations.
* output: CSV file.
"""

import pandas as pd
from scipy.stats import pearsonr, spearmanr


df = snakemake.params['params_df'].copy()
df['pearson'] = 0.
df['spearman'] = 0.


for i, row in df.iterrows():
	df2 = pd.read_csv(snakemake.input[i])

	mask = ~pd.isna(df2['ani'])
	x = df2['gambit'][mask]
	y = df2['ani'][mask]

	df.loc[i, 'pearson'] = pearsonr(x, y)[0]
	df.loc[i, 'spearman'] = spearmanr(x, y)[0]

df.to_csv(snakemake.output[0], index=False)
