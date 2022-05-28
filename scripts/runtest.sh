#!/usr/bin/env bash 

# Run SnakeMake in test mode using $REPO_ROOT/test-workdir/ as the working directory

set -euo pipefail


# https://stackoverflow.com/a/246128/1775059
SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]:-$0}"; )" &> /dev/null && pwd 2> /dev/null; )";
REPO_DIR="$(dirname $(realpath "$SCRIPT_DIR"))"
SNAKEFILE="$REPO_DIR/workflow/Snakefile"
WORKDIR="$REPO_DIR/test-workdir"

set -x
snakemake -s "$SNAKEFILE" -d "$WORKDIR" "$@" --config test=1
