"""
Tools for running benchmarks.
"""

from pathlib import Path
import subprocess as sp
import typing as t
import shlex

from attrs import define, field


def pass_none(f):
	def wrapper(x):
		return None if x is None else f(x)

	return wrapper


@define
class BenchmarkEnv:
	"""Global settings for all benchmarks."""
	ncores: int = field()


@define
class Command:
	args: t.List[str] = field(converter=lambda l: list(map(str, l)))

	def __str__(self):
		return shlex.join(self.args)


@define
class TimeResult:
	real: float
	user: float
	sys: float

	def __str__(self):
		return f'({self.real:.2f} real, {self.user:.2f} user, {self.sys:.2f} sys)'


@define
class BenchmarkResult:
	command: Command
	time: TimeResult


def run(cmd, stdout=sp.DEVNULL, stderr=sp.DEVNULL, **kw) -> sp.CompletedProcess:
	"""Same as subprocess.run() but with different defaults.

	Also accepts path-like for stdout and stderr arguments.
	"""
	if isinstance(cmd, Command):
		cmd = cmd.args

	kw.setdefault('check', True)

	if isinstance(stdout, (str, Path)):
		stdout = open(stdout, 'wb')
	if isinstance(stderr, (str, Path)):
		stderr = open(stderr, 'wb')

	return sp.run(cmd, stdout=stdout, stderr=stderr, **kw)


def benchmark_command(command: Command, workdir, result_file, **kw):
	workdir = Path(workdir)

	args = make_time_cmd(command.args, result_file)
	run(args, cwd=str(workdir), **kw)

	time = parse_time_output(workdir / result_file)
	return BenchmarkResult(
		command=command,
		time=time,
	)


def make_time_cmd(command, result_file):
	"""Make arguments to linux's "time" command."""
	return ['time', '-f', '%e %U %S', '-o', str(result_file), *command]

def parse_time_output(result_file):
	"""Parse output of "time" command."""
	with open(result_file) as f:
		line = f.read()

	parts = line.strip().split()
	assert len(parts) == 3
	return TimeResult(*map(float, parts))


def gambit_signatures_command(list_file, output, params, ncores, relative_to=None):
	args = [
		'gambit', 'signatures', 'create',
		'-c', ncores,
		'-k', params['k'],
		'-p', params['prefix'],
		'-l', list_file,
		'-o', output,
		'--no-progress',
	]
	return Command(args)


def gambit_dist_command(query_sigs, ref_sigs, output, ncores, relative_to=None):
	args = [
		'gambit', 'dist',
		'-c', ncores,
		'--qs', query_sigs,
		'--rs', ref_sigs,
		'-o', output,
		'--no-progress',
	]
	return Command(args)


def fastani_command(query_list, ref_list, output, params, ncores, relative_to=None):
	args = [
		'fastANI',
		'-k', params['k'],
		'--fragLen', params['fraglen'],
		'--threads', ncores,
		'--ql', query_list,
		'--rl', ref_list,
		'-o', output,
	]
	return Command(args)


def mash_sketch_command(files, output, params, relative_to=None):
	# Mash sketching does not support parallelization
	args = [
		'mash', 'sketch',
		'-k', params['k'],
		'-s', params['sketch_size'],
		'-o', output,
		*files
	]
	return Command(args)


def mash_dist_command(query_sketch, ref_sketch, params, ncores, relative_to=None):
	# Output is recorded in stdout
	args = [
		'mash', 'dist',
		'-p', ncores,
		query_sketch,
		ref_sketch,
	]
	return Command(args)
