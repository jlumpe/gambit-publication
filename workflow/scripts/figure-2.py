"""Generate figures 2a and 2b, as well as related supplemental figure 1.

Expected Snakemake variables:

* input: Output of gambit_ani_correlation rule, for appropriate "paramspace" wildcard.
* params
  * paramspace: One of "prefix_length", "prefix_sequence", or "full". This corresponds to which figure
  we're creating.
* output: Figure PNG.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from gambit_pub.utils import genome_set_label


### Plot style ###

plt.style.use('gambit')

plt.rcParams.update({
	'axes.grid': True,
	'axes.grid.axis': 'y',
})

DEFAULT_AXIS_COLOR = 'blue'
DEFAULT_K_LINE_STYLE = dict(
	color=DEFAULT_AXIS_COLOR,
	linestyle='dotted',
)

relplot_kws = dict(
	kind='line',
	height=3,
	dashes=False,
)

FACET_KWS = dict(
	ylim=(.9, 1),
	despine=False,
)


### Load params ###

df = pd.read_csv(snakemake.input[0], index_col=[0, 1, 2, 3])
genome_sets, k_vals, prefix_lens, prefix_versions = df.index.levels

df['genomeset_label'] = list(map(genome_set_label, df.reset_index()['genomeset']))
df['spearman_abs'] = df['spearman'].abs()

# Map (length, version) -> prefix string
prefix_map = df.reset_index().set_index(['prefix_len', 'prefix_version'])['prefix'].drop_duplicates().to_dict()

# Baseline correlation values, corresponding to default choice of parameters
dflt_k = snakemake.config['gambit']['k']
dflt_prefix = snakemake.config['gambit']['prefix']
dflt_plen = len(dflt_prefix)
dflt_pver = 0
assert dflt_prefix == prefix_map[(dflt_plen, dflt_pver)]

### Set plot parameters based on parameter space ###

paramspace = snakemake.params['paramspace']

if paramspace == 'full':
	relplot_kws.update(row='prefix_len', col='prefix_version')
	axis_key_to_prefix = prefix_map
	default_axis_key = [dflt_plen, dflt_pver]

elif paramspace == 'prefix_length':
	relplot_kws['col'] = 'prefix_len'
	axis_key_to_prefix = {
		plen: prefix_map[plen, dflt_pver]
		for plen in prefix_lens
	}
	default_axis_key = dflt_plen

elif paramspace == 'prefix_sequence':
	relplot_kws['col'] = 'prefix_version'
	axis_key_to_prefix = {
		pver: prefix_map[dflt_plen, pver]
		for pver in prefix_versions
	}
	default_axis_key = dflt_pver

else:
	raise ValueError('Incorrect value for "paramspace" job parameter.')


### Plot ###

fg = sns.relplot(
	data=df,
	x='k',
	y='spearman_abs',
	hue='genomeset_label',
	style='genomeset_label',
	markers=['o' for _ in genome_sets],
	facet_kws=FACET_KWS,
	col_wrap=4 if 'row' not in relplot_kws else None,
	**relplot_kws,
)

fg.set_xlabels('$k$')
fg.set_ylabels(r'$|\rho|$', rotation='horizontal', labelpad=20,)
fg.set(xticks=k_vals)
fg.legend.set_title('Genome Set')

# Remove ticks for inner subplots
for ax in set(fg._not_bottom_axes):
	ax.tick_params(bottom=False)
for ax in set(fg._not_left_axes):
	ax.tick_params(left=False)

# Set axis titles to prefix sequences
for key, ax in fg.axes_dict.items():
	ax.set_title(axis_key_to_prefix[key])

# Highlight default axis
dflt_ax = fg.axes_dict[default_axis_key]
for spine in dflt_ax.spines.values():
	spine.set_color(DEFAULT_AXIS_COLOR)
dflt_ax.axvline(dflt_k, **DEFAULT_K_LINE_STYLE)

fg.tight_layout()


### Save ###

fg.figure.savefig(snakemake.output[0])
