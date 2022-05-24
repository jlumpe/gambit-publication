"""Create figure 5.

Expected Snakemake variables:

* input
	db_genomes: Database genomes file.
	db_signatures: Database signatures file.
* params
	conf:
		in_taxon: ID of "in" taxon. Must be divided into subgroups.
		out_taxa: IDs of "out" taxa.
		bin_width: Override default histogram bin width.
* output: Figure png.
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from gambit.db import ReferenceDatabase, Taxon
from gambit.metric import jaccarddist_matrix, jaccarddist_pairwise
from gambit.util.misc import zip_strict


conf = snakemake.params['conf']


### Plot style ###

plt.style.use('gambit')

plt.rcParams.update({
	'axes.grid.axis': 'x',
})

palette = plt.rcParams['axes.prop_cycle'].by_key()['color']

SUBPLOT_SIZE =  (12, 3)
BINWIDTH = conf.get('bin_width', 0.01)

INTRA_COLOR = palette[0]
INTER_SG_COLOR = palette[1]
OUT_PALETTE = palette[2:]

HISTOGRAM_STYLE = dict(
	stat='percent',
	element='step',
	alpha=.1,
)
THRESHOLD_LINE_STYLE = dict(lw=2, linestyle='dashed')


### Setup ###

refdb = ReferenceDatabase.load(snakemake.input['db_genomes'], snakemake.input['db_signatures'])

in_taxon = refdb.session.query(Taxon).get(conf['in_taxon'])
out_taxa = [refdb.session.query(Taxon).get(tid) for tid in conf['out_taxa']]
n_out = len(out_taxa)

subgroups = list(in_taxon.children)
n_subgroups = len(subgroups)

assert len(subgroups) > 0
assert all(sg.isleaf() for sg in subgroups)


### Load signatures ###

genome_tids = np.asarray([genome.taxon_id for genome in refdb.genomes])

# Subgroups
sg_gidxs = [np.flatnonzero(genome_tids == sg.id) for sg in subgroups]
sg_sizes = list(map(len, sg_gidxs))
sg_sigs = [refdb.signatures[gidxs] for gidxs in sg_gidxs]

# "Out" taxa
out_sigs = []

for taxon in out_taxa:
	subtree_tids = [t.id for t in taxon.traverse()]
	gidxs = np.flatnonzero(np.in1d(genome_tids, subtree_tids))
	out_sigs.append(refdb.signatures[gidxs])


### Calculate distances ###

# Within subgroups
subgroup_intra = [jaccarddist_pairwise(sigs, flat=True) for sigs in sg_sigs]

# Between different subgroups
inter_pairs = [(i, j) for i in range(n_subgroups) for j in range(i+1, n_subgroups)]
subgroup_inter = dict()

for i, j in inter_pairs:
	dists = jaccarddist_matrix(sg_sigs[i], sg_sigs[j]).flatten()
	subgroup_inter[i, j] = subgroup_inter[j, i] = dists

subgroup_inter_combined = [
	np.concatenate([subgroup_inter[i, j] for j in range(n_subgroups) if i != j])
	for i in range(n_subgroups)
]

# Within "in" taxon
full_intra = np.concatenate([*subgroup_intra, *(subgroup_inter[p] for p in inter_pairs)])

# Between subgroups and "out" taxa
sg_out_inter = np.empty((n_subgroups, n_out), dtype=object)

for i, sigs1 in enumerate(sg_sigs):
	for j, sigs2 in enumerate(out_sigs):
		sg_out_inter[i, j] = jaccarddist_matrix(sigs1, sigs2).flatten()

# Between "in" and "out" taxa
full_inter = [np.concatenate(col) for col in sg_out_inter.T]

# Some quick sanity checks
assert np.all(full_intra > 0)
for i, sg_taxon in enumerate(subgroups):
	assert all(d.min() > sg_taxon.distance_threshold for d in sg_out_inter[i, :])


### Plot ###

nrow = n_subgroups + 1
fig, axes = plt.subplots(
	nrow, 1,
	figsize=(SUBPLOT_SIZE[0], nrow * SUBPLOT_SIZE[1]),
	sharex=True,
)

xmax = max(full_intra.max(), max(dists.max() for dists in full_inter))
out_colors = OUT_PALETTE[:n_out]

def histogram(ax, x, color, **kw):
	return sns.histplot(
		ax=ax,
		x=x,
		color=color,
		binwidth=BINWIDTH,
		binrange=(0, xmax),
		**HISTOGRAM_STYLE,
		**kw,
	)


# Top subplot - all subgroups combined
histogram(axes[0], full_intra, INTRA_COLOR, label='Intra-species/subgroup')
for taxon, dists, color in zip_strict(out_taxa, full_inter, out_colors):
	histogram(axes[0], dists, color, label=f'Inter-species ({taxon.name})')

axes[0].set_title(in_taxon.name)


# Remaining subplots - individual subgroups
for i, ax in enumerate(axes[1:]):
	ax.set_title(f'Subgroup {i+1} (n={sg_sizes[i]})')

	# Subgroup intra/inter
	histogram(ax, subgroup_intra[i], INTRA_COLOR)
	histogram(ax, subgroup_inter_combined[i], INTER_SG_COLOR, label='Inter-subgroup' if i == 0 else None)

	# Subgroup inter-species
	for dists, color in zip_strict(sg_out_inter[i], out_colors):
		histogram(ax, dists, color)

	# Threshold
	ax.axvline(
		subgroups[i].distance_threshold,
		color=INTRA_COLOR,
		label='Subgroup threshold' if i == 0 else None,
		**THRESHOLD_LINE_STYLE,
	)


# Finish up
axes[-1].set_xlabel('GAMBIT Distance')
plt.tight_layout()
fig.legend(
	loc='upper center',
	ncol=3 + n_out,
	bbox_to_anchor=(.5, 0),
)


### Save ###

fig.savefig(snakemake.output[0], bbox_inches='tight')
