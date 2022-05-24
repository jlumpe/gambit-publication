"""Create one of the figure 4 subplots.

Expected Snakemake variables:

* input
	db_genomes: Database genomes file.
	db_signatures: Database signatures file.
* params
	conf:
		in_taxon: ID of "in" taxon.
		out_taxa: IDs of "out" taxa.
* output: Figure png.
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from gambit.db import ReferenceDatabase, Taxon
from gambit.metric import jaccarddist_matrix, jaccarddist_pairwise


### Plot style ###

plt.style.use('gambit')

palette = plt.rcParams['axes.prop_cycle'].by_key()['color']

FIGSIZE =  (12, 4)
BINWIDTH = .01


### Setup ###

refdb = ReferenceDatabase.load(snakemake.input['db_genomes'], snakemake.input['db_signatures'])

conf = snakemake.params['conf']
in_taxon = refdb.session.query(Taxon).get(conf['in_taxon'])
out_taxa = [refdb.session.query(Taxon).get(tid) for tid in conf['out_taxa']]

in_thresh = in_taxon.distance_threshold
assert in_thresh is not None


### Load signatures ###

genome_tids = np.asarray([genome.taxon_id for genome in refdb.genomes])

in_subtree_tids = [taxon.id for taxon in in_taxon.traverse()]
in_gidxs = np.flatnonzero(np.in1d(genome_tids, in_subtree_tids))
in_ngenomes = len(in_gidxs)
in_sigs = refdb.signatures[in_gidxs]

out_sigs = []

for taxon in out_taxa:
	subtree_tids = [t.id for t in taxon.traverse()]
	gidxs = np.flatnonzero(np.in1d(genome_tids, subtree_tids))
	out_sigs.append(refdb.signatures[gidxs])


### Calculate distances ###

intra_dists = jaccarddist_pairwise(in_sigs)
inter_dists = [jaccarddist_matrix(in_sigs, os) for os in out_sigs]

intra_dists_flat = intra_dists[np.triu_indices(in_ngenomes, k=1)]

# Mask diagonal values before finding minima
mask = np.identity(in_ngenomes, dtype=bool)
intra_masked = np.ma.masked_array(intra_dists, mask)

# argmin_intra = np.argmin(intra_masked, axis=0)
min_intra = np.min(intra_masked, axis=1)

min_inter = [np.min(dmat, axis=1) for dmat in inter_dists]


# Just take this time to check that our database isn't screwed up
for dmat in inter_dists:
	assert np.all(dmat > in_thresh)

for mi in min_inter:
	assert np.all(min_intra < mi)


### Plot ###

fig = plt.figure(figsize=FIGSIZE)
ax = plt.gca()
ax.yaxis.set_visible(False)

values = [intra_dists_flat, min_intra, *min_inter]
labels = ['All intra-species', 'Minimum']
labels.extend([f'Minimum inter-species ({taxon.name})' for taxon in out_taxa])
xmax = max(v.max() for v in values)

# Histograms
for x, label, color in zip(values, labels, palette):
	sns.histplot(
		x=x,
		element='step',
		binrange=(0, xmax),
		binwidth=BINWIDTH,
		alpha=.3,
		color=color,
		stat='density',
		label=label,
	)

# In taxon distance threshold
ax.axvline(
	in_thresh,
	lw=2,
	linestyle='dashed',
	color=palette[0],
	label='Threshold',
)


# Finish up
ax.set_title(in_taxon.name)
ax.set_xlabel('GAMBIT Distance')
plt.tight_layout()
ax.legend(
	ncol=len(labels) + 1,
	# loc='lower center',
	loc='upper right',
)

# Adjust y limits to give some space for the legend
ax.set_ylim(0, ax.get_ylim()[1] * 1.1)


### Save ###

fig.savefig(snakemake.output[0])
