"""
Common code for figure 2 and supplemental figure 1 which explore the main GAMBIT parameter space.
"""

from types import SimpleNamespace

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from gambit_pub.utils import genome_set_label


STYLE_PARAMS = {
	'axes.grid': True,
	'axes.grid.axis': 'y',
}

DEFAULT_AXIS_COLOR = 'blue'

RELPLOT_KWS = dict(
	kind='line',
	height=3,
	dashes=False,
)

FACET_KWS = dict(
	ylim=(.9, 1),
	despine=False,
)


def make_params_df(k, prefix_len, base_prefixes, genome_sets):
	k = np.atleast_1d(k)
	prefix_len = np.atleast_1d(prefix_len)
	base_prefixes = np.atleast_1d(base_prefixes)
	genome_sets = np.atleast_1d(genome_sets)

	index = pd.MultiIndex.from_product(
		[genome_sets, k, prefix_len, range(len(base_prefixes))],
		names=['genomeset', 'k', 'prefix_len', 'prefix_version'],
	)
	df = index.to_frame(False)

	df['base_prefix'] = base_prefixes[df['prefix_version']]
	df['prefix'] = [row.base_prefix[:row.prefix_len] for _, row in df.iterrows()]

	return df


def get_param_data(stats_file, config):
	"""Calculate variables relevant to parameter space to be explored.

	Parameters
	----------
	params_file
		CSV file created by the `gambit-vs-ani-correlation.py` script.
	config
		Snakemake config dict.

	Returns
	-------
	types.SimpleNamespace
		Namespace object containing data.
	"""
	data = SimpleNamespace()

	# Load parameters dataframe
	data.df = pd.read_csv(stats_file, index_col=[0, 1, 2, 3])
	data.genome_sets, data.k_vals, data.prefix_lens, data.prefix_versions = data.df.index.levels

	data.df['genomeset_label'] = list(map(genome_set_label, data.df.index.get_level_values('genomeset')))
	data.df['spearman_abs'] = data.df['spearman'].abs()

	# Map (length, version) -> prefix string
	data.prefix_map = data.df.reset_index().set_index(['prefix_len', 'prefix_version'])['prefix'].drop_duplicates().to_dict()

	# Baseline correlation values, corresponding to default choice of parameters
	data.dflt_k = config['gambit']['k']
	data.dflt_prefix = config['gambit']['prefix']
	data.dflt_plen = len(data.dflt_prefix)
	data.dflt_pver = 0
	assert data.dflt_prefix == data.prefix_map[(data.dflt_plen, data.dflt_pver)]

	return data


def set_style():
	"""Set Matplotlib style variables."""
	plt.style.use('gambit')
	plt.rcParams.update(STYLE_PARAMS)


def spearman_vs_k(paramdata, **kw):
	"""
	Plot (absolute) Spearman correlation vs. k colored by genome set, possibly split over subplots.
	"""

	facet_kws = {**FACET_KWS, **kw.pop('facet_kws', {})}
	kw = {**RELPLOT_KWS, **kw}

	fg = sns.relplot(
		data=paramdata.df,
		x='k',
		y='spearman_abs',
		hue='genomeset_label',
		style='genomeset_label',
		markers=['o' for _ in paramdata.genome_sets],
		facet_kws=facet_kws,
		**kw,
	)

	fg.set_xlabels('$k$')
	fg.set_ylabels('$\\rho$', rotation='horizontal')
	fg.set(xticks=paramdata.k_vals)
	fg.legend.set_title('Genome Set')

	# Remove ticks for inner subplots
	for ax in set(fg._not_bottom_axes):
		ax.tick_params(bottom=False)
	for ax in set(fg._not_left_axes):
		ax.tick_params(left=False)

	return fg


def highlight_default_axis(ax):
	for spine in ax.spines.values():
		spine.set_color(DEFAULT_AXIS_COLOR)
