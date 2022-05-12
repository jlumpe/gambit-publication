"""Generate figure 1.

Expected Snakemake variables:

* input: gambit_vs_ani rule output for each genome set.
* params
  * genome_sets: Genome set IDs.
* output:
  * figure: Figure PNG.
  * stats: Statistics CSV file.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr
from matplotlib.ticker import PercentFormatter
from gambit.util.misc import zip_strict


# ## Setup

genome_sets = snakemake.params['genome_sets']
ngs = len(genome_sets)

mpl_config = snakemake.config['matplotlib']


### Plot style ###

plt.style.use(mpl_config['style'])

plt.rcParams.update({
	'axes.grid': True,
})

SUBPLOT_SIZE = 5
NBINS = 100


### Load data and do calculations ###

plot_data = []
stats_rows = []

for gset, file in zip_strict(genome_sets, snakemake.input):
	data = dict()

	df = pd.read_csv(file, index_col=[0, 1])
	df_filtered = df.dropna()
	corr = spearmanr(df_filtered['ani'], df_filtered['gambit']).correlation

	stats_rows.append((gset, df.shape[0], df_filtered.shape[0], corr))
	plot_data.append(df_filtered)

stats_df = pd.DataFrame(stats_rows, columns=['genome_set', 'npairs', 'ani_reported', 'spearmanr'])


### Plot ###

fig, axes = plt.subplots(
	1, ngs,
	figsize=(ngs * SUBPLOT_SIZE, SUBPLOT_SIZE),
	sharex=True,
	sharey=True,
)

for ax, (_, row), df in zip_strict(axes, stats_df.iterrows(), plot_data):
	ax.set_title(row.genome_set)
	ax.set_xlabel('ANI')

	# Spearman correlation
	ax.text(
		.05, .05,
		f'$\\rho = {row.spearmanr:.3f}$',
		ha='left',
		va='bottom',
		fontsize=16,
		transform=ax.transAxes,
	)

	# ax.scatter(df['ani'], df['gambit'], s=1)
	sns.histplot(
		data=df,
		x='ani',
		y='gambit',
		ax=ax,
		bins=NBINS,
		binrange=(
			(df['ani'].min(), 100),
			(0, 1),
		),
	)


### Finish axes ####

axes[0].set_ylabel('GAMBIT Distance')
axes[0].set_xlim(None, 100)
axes[0].set_ylim(0, 1)
axes[0].xaxis.set_major_formatter(PercentFormatter())

fig.tight_layout()


### Save ###

fig.savefig(snakemake.output['figure'], **mpl_config['savefig_args'])

stats_df.to_csv(snakemake.output['stats'], index=False)
