"""
Tools for running benchmarks.
"""

import sys
from pathlib import Path
import subprocess as sp
import typing as t
import shlex
import random

from attrs import define, field

from .utils import truncate_str


def pass_none(f):
	def wrapper(x):
		return None if x is None else f(x)

	return wrapper


@define
class BenchmarkEnv:
	"""Global settings for all benchmarks."""
	ncores: int = field()


@define
class TimeResult:
	real: float
	user: float
	sys: float

	def __str__(self):
		return f'({self.real:.2f} real, {self.user:.2f} user, {self.sys:.2f} sys)'


def run(cmd, stdout=sp.DEVNULL, stderr=sp.DEVNULL, **kw) -> sp.CompletedProcess:
	"""Same as subprocess.run() but with different defaults.

	Also accepts path-like for stdout and stderr arguments.
	"""
	kw.setdefault('check', True)

	if isinstance(stdout, (str, Path)):
		stdout = open(stdout, 'wb')
	if isinstance(stderr, (str, Path)):
		stderr = open(stderr, 'wb')

	return sp.run(cmd, stdout=stdout, stderr=stderr, **kw)


def time_command(cmd, workdir, result_file, **kw) -> TimeResult:
	workdir = Path(workdir)

	args = make_time_cmd(cmd, result_file)
	run(args, cwd=str(workdir), **kw)

	return parse_time_output(workdir / result_file)


def make_time_cmd(cmd, result_file):
	"""Make arguments to linux's "time" command."""
	return ['time', '-f', '%e %U %S', '-o', str(result_file), *cmd]

def parse_time_output(result_file):
	"""Parse output of "time" command."""
	with open(result_file) as f:
		line = f.read()

	parts = line.strip().split()
	assert len(parts) == 3
	return TimeResult(*map(float, parts))


@define
class BenchmarkResult:
	# command: Command
	key: t.Tuple
	command: t.List[str]
	workdir: Path
	round: int
	execution_order: int
	time: TimeResult


class BenchmarkRunner:
	"""Runs a set of commands multiple times and records the execution times."""

	def __init__(self, workdir='.', log_file=sys.stderr):
		self.workdir = Path(workdir)
		self.log_file = log_file
		self.commands = dict()
		self.truncate_print = 400

	def add_command(self, key, cmd: t.List[str]):
		if not isinstance(key, tuple):
			key = (key,)
		if key in self.commands:
			raise KeyError(f'Command with key {key!r} already exists')
		self.commands[key] = list(map(str, cmd))

	def _key_to_str(self, key):
		return '-'.join(map(str, key))

	def make_cmd_dir(self, key, r):
		d = self.workdir / Path(f'{self._key_to_str(key)}-{r}')
		d.mkdir()
		return d

	def log(self, *args, leading=None, trailing=None, **kw):
		"""Display a log message."""
		if leading is not None:
			print('\n' * leading, file=self.log_file, end='')
		print(*args, file=sys.stderr, **kw)
		if trailing is not None:
			print('\n' * trailing, file=self.log_file, end='')

	def log_command(self, cmd: t.List[str], **kw):
		"""Log command before executing it."""
		s = shlex.join(cmd)
		s = truncate_str(s, self.truncate_print)
		self.log(f'+ {s}', **kw)

	def run(self,
			replicates: int,
			dry_run: bool = False,
	        shuffle: bool = True,
			) -> t.List[BenchmarkResult]:

		self.workdir.mkdir(exist_ok=True)

		results = []

		# Each round
		for r in range(replicates):
			self.log(f'*** Starting round {r + 1} of {replicates} ***', leading=2)

			# Shuffle list of commands
			keys_list = list(self.commands)
			if shuffle:
				random.seed(r)
				random.shuffle(keys_list)

			# Each command
			for i, key in enumerate(keys_list):
				keystr = self._key_to_str(key)
				cmd = self.commands[key]

				self.log(f'Running {keystr} ({i + 1} of {len(self.commands)})', leading=1)
				self.log_command(cmd)

				if dry_run:
					continue

				workdir = self.make_cmd_dir(key, r)

				# Try running
				try:
					time = time_command(
						cmd, workdir, 'time-out',
						stdout=workdir / 'stdout',
						stderr=workdir / 'stderr',
					)

				except sp.CalledProcessError as e:
					# Print stdout and stderr and re-raise
					with open(workdir / 'stdout') as f:
						self.log(f.read())
					with open(workdir / 'stderr') as f:
						self.log(f.read())
					raise

				self.log(f'Completed in {time}')

				results.append(BenchmarkResult(key, cmd, workdir, r, i, time))

		return results
