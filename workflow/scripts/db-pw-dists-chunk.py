"""Calculate distance matrix between two of reference genome chunks.

Expected Snakemake variables:

* input
	db_signatures: Database signatures file.
	chunks_dir: Directory containing genome tables for each chunk.
* wildcards
	chunk1: Index of 1st chunk.
	chunk2: Index of 2nd chunk (>= first).
* params
	* show_progress: Whether to display the progress bar
* output
"""

from pathlib import Path

import pandas as pd
import h5py as h5

from gambit.sigs import load_signatures
from gambit.metric import jaccarddist_matrix


### Setup ###

chunk1 = int(snakemake.wildcards['chunk1'])
chunk2 = int(snakemake.wildcards['chunk2'])
# Make sure we're keeping the correct (arbitrary) order for chunk pairs
assert chunk1 <= chunk2

chunks_dir = Path(snakemake.input['chunks_dir'])
chunk1_genomes = pd.read_csv(chunks_dir / f'{chunk1}.csv')
chunk2_genomes = pd.read_csv(chunks_dir / f'{chunk2}.csv')

show_progress = snakemake.params.get('show_progress', False)


### Calculate distances ###

sigs = load_signatures(snakemake.input['signatures'])

chunk1_sigs = sigs[chunk1_genomes['signatures_index']]
chunk2_sigs = sigs[chunk2_genomes['signatures_index']]

dmat = jaccarddist_matrix(chunk1_sigs, chunk2_sigs, progress=show_progress)


### Save ###

with h5.File(snakemake.output[0], 'w') as f:
	f.create_dataset('dmat', data=dmat)
