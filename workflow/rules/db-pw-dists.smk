"""Rules for generating pairwise distance matrix for GAMBIT database reference genomes."""


# Divide reference genomes into chunks based on taxonomy
checkpoint db_pw_dists_make_chunks:
	input:
		db_signatures=rules.get_gambit_db.output['signatures'],
		db_genomes=rules.get_gambit_db.output['genomes'],
	params:
		max_chunk_size=2500,
	output:
		summary_table='intermediate-data/db-pw-dists/chunks.csv',
		taxa_table='intermediate-data/db-pw-dists/chunk-taxa.csv',
		chunks_dir=directory('intermediate-data/db-pw-dists/chunks/'),
	script: '../scripts/db-pw-dists-make-chunks.py'


# Create pairwise distance matrix for two chunks
rule db_pw_dists_chunk:
	input:
		signatures=rules.get_gambit_db.output['signatures'],
		chunks_dir=rules.db_pw_dists_make_chunks.output['chunks_dir'],
	output: 'intermediate-data/db-pw-dists/dists/{chunk1}-{chunk2}.h5'
	wildcard_constraints:
		chunk1='\d+',
		chunk2='\d+',
	params:
		show_progress=config['show_progress'],  # Display progress bar
	threads: workflow.cores
	script: '../scripts/db-pw-dists-chunk.py'


def get_db_pw_dists_genome_chunks(wildcards=None):
	cp = checkpoints.db_pw_dists_make_chunks.get()
	chunks_df = pd.read_csv(cp.output['summary_table'], index_col=0)
	nchunks = chunks_df.shape[0]
	assert np.array_equal(chunks_df.index, range(nchunks))

	chunks_dir = Path(cp.output['chunks_dir'])
	return [chunks_dir / f'{i}.csv' for i in range(nchunks)]


def get_db_pw_dists_dmat_chunks(wildcards=None):
	chunks_file = checkpoints.db_pw_dists_make_chunks.get().output['summary_table']
	chunks_df = pd.read_csv(chunks_file, index_col=0)
	nchunks = chunks_df.shape[0]
	assert np.array_equal(chunks_df.index, range(nchunks))

	files = []

	for i in range(nchunks):
		for j in range(i, nchunks):
			files.extend(expand(rules.db_pw_dists_chunk.output, chunk1=i, chunk2=j))

	return files


# Minimum and maximum distances between leaf taxa
rule db_taxa_pw_dist_ranges:
	input:
		db_genomes=rules.get_gambit_db.output['genomes'],
		genome_chunks=get_db_pw_dists_genome_chunks,
		dmat_chunks=get_db_pw_dists_dmat_chunks,
	output:
		min_dists='intermediate-data/db-pw-dists/leaf-min-dists.csv',
		max_dists='intermediate-data/db-pw-dists/leaf-max-dists.csv',
	script: '../scripts/db-taxa-pw-dist-ranges.py'
