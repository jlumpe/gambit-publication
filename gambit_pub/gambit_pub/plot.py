"""Plotting tools."""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
import seaborn as sns


def get_subplots_bbox(axes: np.ndarray) -> Bbox:
	"""Get bounding box of subplots in figure coordinates."""
	xmin = np.inf
	xmax = -np.inf
	ymin = np.inf
	ymax = -np.inf

	for ax in axes.flat:
		bb = ax.get_position()
		xmin = min(xmin, bb.xmin)
		xmax = max(xmax, bb.xmax)
		ymin = min(ymin, bb.ymin)
		ymax = max(ymax, bb.ymax)

	return Bbox.from_extents(xmin, ymin, xmax, ymax)


def subplots_xlabel(fig: plt.Figure, axes: np.ndarray, label: str, **kw):
	"""Add single common x label for subplots."""
	bbox = get_subplots_bbox(axes)
	fig.supxlabel(label, x=(bbox.xmin + bbox.xmax) / 2, **kw)

def subplots_ylabel(fig: plt.Figure, axes: np.ndarray, label: str, **kw):
	"""Add single common x label for subplots."""
	bbox = get_subplots_bbox(axes)
	fig.supylabel(label, y=(bbox.ymin + bbox.ymax) / 2, **kw)


def set_fg_xlabel(fg: sns.FacetGrid, label: str, **kw):
	"""Add single common x label for FacetGrid subplots."""
	fg.set_xlabels('')
	subplots_xlabel(fg.figure, fg.axes, label, **kw)

def set_fg_ylabel(fg: sns.FacetGrid, label: str, **kw):
	"""Add single common y label for FacetGrid subplots."""
	fg.set_ylabels('')
	subplots_ylabel(fg.figure, fg.axes, label, **kw)


def remove_fg_inner_ticks(fg: sns.FacetGrid, x=True, y=True):
	"""Hide tick marks on inner shared axes of FacetGrid."""

	if x:
		for ax in set(fg._not_bottom_axes):
			ax.tick_params(bottom=False)
	if y:
		for ax in set(fg._not_left_axes):
			ax.tick_params(left=False)
