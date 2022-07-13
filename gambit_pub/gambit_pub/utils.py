"""Misc utility functions."""

from pathlib import Path
import os
import re

import pandas as pd


def symlink_to_relative(dst, src, is_dir=False):
	"""Create a symlink frm src to dst, where both paths are relative to the same directory.

	This is different than os.symlink and pathlib.Path.symlink_to, in which dst is taken to be
	relative to the directory of src.
	"""
	src = Path(src)
	dst_rel = os.path.relpath(dst, src.parent)
	src.symlink_to(dst_rel, is_dir)


def genome_set_label(gset):
	"""Get a nicely formatted label for a genome set given its wildcard value."""
	m = re.fullmatch(r'set(\d[a-z]?)', gset)
	assert m is not None
	return 'Set ' + m.group(1).upper()


def truncate_str(s, n):
	l = len(s)
	if l <= n:
		return s
	return f'{s[:n]}... ({l - n} characters omitted)'


def getattr_coalesce(x, *attrnames: str):
	"""Get attribute value from object ``x``, returning None if ``x`` is None."""
	for name in attrnames:
		if x is None:
			return None
		x = getattr(x, name)
	return x


def fix_int_cols(df: pd.DataFrame, cols):
	"""Fix DataFrame columns which were automatically converted to float due to null values."""
	if isinstance(cols, str):
		cols = [cols]

	for col in cols:
		df[col] = df[col].astype(pd.Int64Dtype())
