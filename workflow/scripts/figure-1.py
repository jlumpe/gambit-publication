"""Generate figure 1.

Expected Snakemake variables:

* input:
  * pairs: "pairs" output of gambit_vs_ani rule for each genome set.
  * stats: "stats" output of gambit_vs_ani rule for each genome set.
* params
  * genome_sets: Genome set IDs.
* output:
  * figure: Figure PNG.
  * stats: Statistics CSV file.
"""

import json

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
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

for gset, pairs_in, stats_in in zip_strict(genome_sets, snakemake.input['pairs'], snakemake.input['stats']):
	df = pd.read_csv(pairs_in, index_col=[0, 1])
	plot_data.append(df.dropna())

	with open(stats_in) as f:
		stats_rows.append(json.load(f))

plot_df = pd.concat(plot_data, keys=genome_sets, names=['genome_set'])

stats_df = pd.DataFrame(stats_rows, index=pd.Series(genome_sets, name='genome_set'))


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
	rho = stats_df.loc[gset, 'spearman']
	ax.text(
		.05, .05,
		f'$\\rho = {rho:.3f}$',
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
