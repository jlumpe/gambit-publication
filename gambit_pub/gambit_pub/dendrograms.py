import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection


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
