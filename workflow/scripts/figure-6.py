"""Generate figure 6.

Expected Snakemake variables:

* input
    * genomes_csv: CSV file of genome attributes.
    * pw_dists: CSV file of pairwise GAMBIT distances.
* output: Figure PNG
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.collections import LineCollection
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage


### Plot style ###

mpl_config = snakemake.config['matplotlib']

plt.style.use(mpl_config['style'])

PHYLOGROUP_PALETTE = 'Set1'
FIG_HEIGHT = 4
CMAP = 'Purples_r'


### Code ###

def linkage_to_df(link):
    """Convert array returned by scipy...linkage() to DataFrame of node attributes."""
    nleaves = link.shape[0] + 1
    nnodes = 2 * nleaves - 1

    # Fill with leaf values
    df = pd.DataFrame.from_dict(dict(
        left=np.full(nnodes, -1),
        right=np.full(nnodes, -1),
        height=np.full(nnodes, 0.),
        size=np.full(nnodes, 1),
    ))

    # Copy values for internal nodes from linkage array
    for i in range(4):
        df.iloc[:, i].values[nleaves:] = link[:, i]

    return df

def make_dendrogram(link_df, left=0.):
    """Make a DataFrame with coordinates for dendrogram branches."""
    nnodes = link_df.shape[0]
    nleaves = (nnodes + 1) // 2

    nodes = link_df[['left', 'right', 'height', 'size']]
    nodes['center'] = 0
    _make_subtree(nodes, nodes.index[-1], 0)

    leaf_order = np.empty(nleaves, dtype=int)
    for i in range(nleaves):
        c = nodes.loc[i, 'center']
        assert float(c).is_integer()
        leaf_order[int(c)] = i

    nodes['center'] += left

    return dict(
        nodes=nodes,
        leaf_order=leaf_order,
    )

def _make_subtree(nodes, i, ll):
    """Recursively fill in info for nodes, bottom to top. Return (center, right_leaf)."""

    if nodes.loc[i, 'size'] <= 1:
        nodes.loc[i, 'center'] = ll
        return (ll, ll)

    else:
        lc, lr = _make_subtree(nodes, nodes.loc[i, 'left'], ll)
        rc, rr = _make_subtree(nodes, nodes.loc[i, 'right'], lr + 1)

        center = nodes.loc[i, 'center'] = (ll + rr) / 2

        return center, rr

def draw_dendrogram(ax, dg, horizontal=False, colorfunc=None):
    """Draw a dendrogram onto a matplotlib axes object."""
    nodes = dg['nodes']
    nnodes = nodes.shape[0]
    nleaves = (nnodes + 1) // 2
    internal = range(nleaves, nnodes)

    segments = [_node_segment(nodes, i, horizontal) for i in internal]
    colors = [colorfunc(i) for i in internal]

    lc = LineCollection(segments, colors=colors)
    ax.add_collection(lc)
    ax.autoscale()

    if horizontal:
        ax.yaxis.set_visible(False)
    else:
        ax.xaxis.set_visible(False)

def _node_segment(nodes, i, horizontal):
    li, ri, h = nodes.loc[i, ['left', 'right', 'height']]
    lh, lc = nodes.loc[li, ['height', 'center']]
    rh, rc = nodes.loc[ri, ['height', 'center']]

    segment = [(lc, lh), (lc, h), (rc, h), (rc, rh)]
    return [p[::-1] for p in segment] if horizontal else segment


### Load input data ###

genomes_df = pd.read_csv(snakemake.input['genomes_csv'])
ngenomes = genomes_df.shape[0]

dists_df = pd.read_csv(snakemake.input['pw_dists'], index_col=0)
assert np.array_equal(dists_df.index, dists_df.columns)
assert np.array_equal(dists_df.index, [id_ + '.fasta' for id_ in genomes_df['id']])
dmat = dists_df.values


### Clustering ###

link = linkage(squareform(dmat), 'average')
nodes  = linkage_to_df(link)
nnodes = nodes.shape[0]

# Assign phylogroups
nodes.loc[:ngenomes, 'phylogroup'] = genomes_df['phylogroup']
for i in range(ngenomes, nnodes):
    left_pg = nodes.loc[nodes.loc[i, 'left'], 'phylogroup']
    right_pg = nodes.loc[nodes.loc[i, 'right'], 'phylogroup']
    if pd.notnull(left_pg) and left_pg == right_pg:
        nodes.loc[i, 'phylogroup'] = left_pg


### Calculated plot parameters ###

phylogroups = genomes_df['phylogroup'].unique()
phylogroups.sort()

phylo_palette = sns.color_palette(PHYLOGROUP_PALETTE, len(phylogroups))
phylo_colors = dict(zip(phylogroups, phylo_palette))
genome_colors = [phylo_colors[pg] for pg in genomes_df['phylogroup']]

# Width of subplots and horizontal padding between them
ax_pad = [.15, .45, .2]
ax_width = [.5 * FIG_HEIGHT, .25, FIG_HEIGHT, .2]

# Left side of axes in physical coordinates (inches)
ax_left = []
for i in range(4):
    left = 0 if i == 0 else ax_left[i-1] + ax_width[i-1] + ax_pad[i-1]
    ax_left.append(left)

fig_w = sum(ax_pad) + sum(ax_width)

# Subplot rects as proportion of figure size (argument to Figure.add_axes())
ax_rects = [(l / fig_w, 0, w / fig_w, 1) for l, w in zip(ax_left, ax_width)]


### Figure/axes ###

fig = plt.figure(figsize=(fig_w, FIG_HEIGHT))
dg_ax = fig.add_axes(ax_rects[0])
tbl_ax = fig.add_axes(ax_rects[1], sharey=dg_ax)
hm_ax = fig.add_axes(ax_rects[2], sharey=dg_ax)
cbar_ax = fig.add_axes(ax_rects[3])


### Dendrogram ###

dg = make_dendrogram(nodes, left=.5)
draw_dendrogram(
    dg_ax,
    dg,
    horizontal=True,
    colorfunc=lambda i: phylo_colors.get(nodes.loc[i, 'phylogroup'], 'black'),
)

dg_ax.invert_xaxis()
dg_ax.invert_yaxis()
dg_ax.set_xlim(None, 0)
dg_ax.set_xlabel('GAMBIT Distance')

for side in ['left', 'top', 'right']:
    dg_ax.spines[side].set_visible(False)


### Heatmap ###

lo = dg['leaf_order']
hm = hm_ax.pcolor(dmat[np.ix_(lo, lo)], cmap=CMAP)
hm_ax.axis('off')
hm_ax.set_aspect(1, share=False, adjustable='box', anchor='W')


### Table ###

for spine in tbl_ax.spines.values():
    spine.set_visible(False)

for i, gi in enumerate(lo):
    pg = genomes_df.loc[gi, 'phylogroup']
    tbl_ax.text(0, i + .5, pg, ha='left', va='center', color=phylo_colors[pg])
    tbl_ax.text(1, i + .5, genomes_df.loc[gi, 'mlst'], ha='left', va='center')

tbl_ax.set_xticks([0, 1])
tbl_ax.set_xticklabels(['Phylogroup', 'MLST'], rotation=-45, ha='left', va='top')
tbl_ax.set_yticks([])
tbl_ax.tick_params(length=0)
tbl_ax.set_xlim(.0, 1.0)


### Color bar ###

plt.colorbar(
    hm,
    cax=cbar_ax,
    label='GAMBIT Distance',
)


### Save ###

fig.savefig(snakemake.output[0], bbox_inches='tight', **mpl_config['savefig_args'])
