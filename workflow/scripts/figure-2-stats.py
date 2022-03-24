"""Generate table of GAMBIT distance/ANI correlations that will be displayed in figure 2."""

import pandas as pd
from scipy.stats import pearsonr, spearmanr


df = snakemake.params['params_df'].copy()
df['pearson'] = 0.
df['spearman'] = 0.

for i, row in df.iterrows():
	df2 = pd.read_csv(snakemake.input[row.input_key])

	mask = ~pd.isna(df2['ani'])
	x = df2['gambit'][mask]
	y = df2['ani'][mask]

	df.loc[i, 'pearson'] = pearsonr(x, y)[0]
	df.loc[i, 'spearman'] = spearmanr(x, y)[0]

df.to_csv(snakemake.output['stats'], index=False)
