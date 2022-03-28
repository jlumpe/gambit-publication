"""Misc utility functions."""

from pathlib import Path
import os


def symlink_to_relative(dst, src, is_dir=False):
    """Create a symlink frm src to dst, where both paths are relative to the same directory.

    This is different than os.symlink and pathlib.Path.symlink_to, in which dst is taken to be
    relative to the directory of src.
    """
    src = Path(src)
    dst_rel = os.path.relpath(dst, src.parent)
    src.symlink_to(dst_rel, is_dir)
