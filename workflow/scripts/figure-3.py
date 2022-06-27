"""Create figure 3.

Expected Snakemake variables:

* input: Directories of fastq_kmers rule output.
* params
	* genomes: Genome IDs.
	* min_phred: Minimum k-mer PHRED score to filter by.
* output: Figure png.
"""

from pathlib import Path
import json
from string import ascii_uppercase

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import gambit_pub.plot as gpp


### Plot style ###

plt.style.use('gambit')

plt.rcParams.update({
	'axes.grid.axis': 'x',
	'axes.grid.which': 'both',
})

palette = plt.rcParams['axes.prop_cycle'].by_key()['color']

HISTOGRAM_STYLE = dict(
	stat='count',
	element='step',
	alpha=.3,
)

SUBPLOT_WIDTH = 8
SUBPLOT_HEIGHT = 2

COVERAGE_LINE_STYLE = dict(lw=1, linestyle='dashed', color='black')
MAJOR_GRID_STYLE = dict(color='#aaaaaa')
MINOR_GRID_STYLE = dict(color='#dddddd')


### Load data ###

ngenomes = len(snakemake.input)
genomes = snakemake.params['genomes']
assert len(genomes) == ngenomes
genome_labels = list(ascii_uppercase[:ngenomes])

min_phred_col = str(snakemake.params['min_phred'])

data_parts = []
coverages = []

for in_dir in map(Path, snakemake.input):
	df = pd.read_csv(in_dir / 'kmer-counts.csv')
	df = df[['in_fasta', min_phred_col]]
	df.columns = ['in_fasta', 'count']
	data_parts.append(df)

	with open(in_dir / 'stats.json') as f:
		stats = json.load(f)

	coverages.append(stats['estimated_coverage'])

data = pd.concat(data_parts, keys=genome_labels, names=['genome_label']).reset_index()
data['K-mer in assembly'] = np.where(data['in_fasta'], 'Yes', 'No')
# import IPython; IPython.embed()

coverages = dict(zip(genome_labels, coverages))


### Histogram bins ###

# max_count = data['count'].max()
max_count = 1000   # Ignore the few outliers outside this range

bins_lin = range(1, 10)
bin_width_log = .1
bins_log = np.arange(1, np.log10(max_count + 1) + bin_width_log, bin_width_log)
bins = np.concatenate([bins_lin, 10 ** bins_log]) - .5


### Plot ###

fg = sns.displot(
	data=data,
	x='count',
	row='genome_label',
	hue='K-mer in assembly',
	bins=bins,
	row_order=genome_labels,
	common_norm=False,
	height=SUBPLOT_HEIGHT,
	aspect=SUBPLOT_WIDTH / SUBPLOT_HEIGHT,
	facet_kws=dict(
		sharey=False,
		margin_titles=True,
		despine=False,
	),
	**HISTOGRAM_STYLE,
)

fg.set_titles(row_template='{row_name}')
gpp.set_fg_xlabel(fg, 'Copy Number')
gpp.set_fg_ylabel(fg, 'Number of k-mers')

# Estimated coverage
for genome, ax in fg.axes_dict.items():
	ax.axvline(
		coverages[genome],
		label='Estimated Coverage',
		**COVERAGE_LINE_STYLE,
	)

# Symlog scale to discrete bins in the 1-10 range
bottom = fg.axes[-1, 0]
bottom.set_xscale('symlog', linthresh=10.5)

# Manually adjust X ticks for the weird scale
xmax = bottom.get_xlim()[1]
maj_ticks_log = np.arange(int(np.ceil(np.log10(xmax))))
maj_ticklabels = ['$1$'] + [f'$10^{n}$' for n in maj_ticks_log[1:]]
min_ticks = [x * 10**n for n in maj_ticks_log for x in range(1, 10)]

for ax in fg.axes.flat:
	# Tick locations/labels
	ax.set_xticks(10 ** maj_ticks_log)
	ax.set_xticklabels(maj_ticklabels)
	ax.set_xticks(min_ticks, minor=True)

	# Grid style
	ax.set_axisbelow(True)
	ax.grid(axis='x', which='major', **MAJOR_GRID_STYLE)
	ax.grid(axis='x', which='minor', **MINOR_GRID_STYLE)

	# Remove ticks on all but bottom subplot
	if ax is not bottom:
		ax.tick_params('x', which='both', bottom=False)


bottom.set_xlim([-0.5, max_count * 1.2])
fg.tight_layout()


### Save ###

plt.savefig(snakemake.output[0])#, bbox_inches='tight')
