#!/usr/bin/env python3

import argparse
from subprocess import run


parser = argparse.ArgumentParser(
	description="Install an IPython kernel spec for this repo's conda environment so that notebooks can be run interactively.",
	formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)

parser.add_argument('--name', default='gambit-publication', help='Name for the kernelspec.')
parser.add_argument('--display-name', default='GAMBIT Publication', help='Human-readable display name.')
parser.add_argument('--no-user', action='store_true', help='Install system-wide instead of for the current user.')
parser.add_argument('-d', '--dry-run', action='store_true', help='Print the installation command instead of running it.')
parser.add_argument('options', nargs='*', help='Additional arguments to pass to the "ipykernel install" command.')


if __name__ == '__main__':
	args = parser.parse_args()

	cmd = ['python', '-m', 'ipykernel', 'install', '--name', args.name, '--display-name', args.display_name]

	if not args.no_user:
		cmd.append('--user')

	cmd.extend(args.options)

	if args.dry_run:
		print(cmd)
	else:
		run(cmd, check=True)
