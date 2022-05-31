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
from gambit_pub.utils import genome_set_label


# ## Setup

genome_sets = snakemake.params['genome_sets']
ngs = len(genome_sets)


### Plot style ###

plt.style.use('gambit')

XRANGE = (75, 100)
YRANGE = (0, 1)

DISPLOT_KW = dict(
	col_wrap=2,
	height=4,
	bins=100,
	binrange=(XRANGE, YRANGE),
	common_norm=False,
	facet_kws=dict(despine=False),
)

CORRELATION_TEXT_KW = dict(
	fontsize=14,
)


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

plot_df = pd.concat(plot_data, keys=genome_sets, names=['genome_set'])

stats_df = pd.DataFrame(stats_rows, columns=['genome_set', 'npairs', 'ani_reported', 'spearmanr'])
stats_df.set_index('genome_set', inplace=True)


### Plot ###
fg = sns.displot(
	data=plot_df.reset_index(),
	col='genome_set',
	x='ani',
	y='gambit',
	**DISPLOT_KW,
)

for gset, ax in fg.axes_dict.items():
	ax.set_title(genome_set_label(gset))

	# Spearman correlation
	rho = stats_df.loc[gset, 'spearmanr']
	ax.text(
		.05, .05,
		f'$\\rho = {-rho:.3f}$',
		ha='left',
		va='bottom',
		transform=ax.transAxes,
		**CORRELATION_TEXT_KW,
	)

	ax.set_xlim(*XRANGE)
	ax.set_ylim(*YRANGE)
	ax.xaxis.set_major_formatter(PercentFormatter(decimals=0))

fg.set_axis_labels('ANI', 'GAMBIT Distance')
fg.tight_layout()


### Save ###

fg.savefig(snakemake.output['figure'])

stats_df.to_csv(snakemake.output['stats'])
