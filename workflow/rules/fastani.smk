"""
Compute all pairwise FastANI comparisons for a given genome set.

Genome set 1 is much larger than the others and FastANI starts to get very slow with this many
genomes. To compensate we divide the query sequences into chunks so that the number of ANI values to
be computed (nqueries * nrefs) is capped to a maximum value (set in config[fastani][chunk_size]).
We run FastANI for each query chunk and then glue the outputs together.
"""


FASTANI_CHUNKS_DIR = 'intermediate-data/fastani/chunks/{genomeset}'


def get_fastani_input(wildcards):
	"""Get input files for the fastani rule (output chunks to combine).

	This depends on the output of the fastani_chunks checkpoint.
	"""
	chunkdir = Path(checkpoints.fastani_query_chunks.get(genomeset=wildcards.genomeset).output[0])
	chunks = [f.stem for f in chunkdir.glob('*-*.txt')]
	return expand(rules.fastani_chunk.output[0], genomeset=wildcards.genomeset, chunk=chunks)


# Divide query genome list into chunks so that the FastANI run time for each is reasonable.
checkpoint fastani_query_chunks:
	input:
		list_file=get_genomes_list_file,
	output: directory(f"{FASTANI_CHUNKS_DIR}/queries")
	run:
		out_dir = Path(output[0])
		out_dir.mkdir()

		# Ready query list file
		with open(input['list_file']) as f:
			lines = list(f)

		# Get chunk size from config (if present), default to single chunk
		chunk_nelems = config['fastani'].get('chunk_size')
		if chunk_nelems is None:
			chunk_nlines = len(lines)
		else:
			nchunks = max(len(lines) ** 2 // chunk_nelems, 1)
			chunk_nlines = len(lines) // nchunks

		# Write chunked list files
		for begin in range(0, len(lines), chunk_nlines):
			chunk = lines[begin:begin + chunk_nlines]
			end = begin + len(chunk)

			with open(out_dir / f'{begin+1:03d}-{end:03d}.txt', 'w') as f:
				f.writelines(chunk)


# Generate single "chunk" of FastANI results, using all genomes for references but a subset of
# genomes for queries.
rule fastani_chunk:
	input:
		fasta=get_genomes_fasta_dir,
		query_chunks=rules.fastani_query_chunks.output[0],
		refs=get_genomes_list_file,
	output: temporary(f'{FASTANI_CHUNKS_DIR}/{{chunk}}.tsv')
	wildcard_constraints:
		chunk="\d+-\d+"
	params:
		k=config['fastani']['k'],
		fraglen=config['fastani']['fraglen'],
	threads: workflow.cores
	shell:
		"""
		 output="$(realpath {output})"
		 queries="$(realpath {input[query_chunks]}/{wildcards.chunk}.txt)"
		 refs="$(realpath {input[refs]})"
		 cd {input[fasta]}
		 fastANI -k {params.k} --fragLen {params.fraglen} -t {threads} \\
			 --ql "$queries" --rl $refs -o "$output"
		 """


# Concatenate chunks to get the complete result file.
rule fastani:
	input: get_fastani_input
	output: protected('intermediate-data/fastani/{genomeset}.tsv')
	shell: 'cat {input} > {output}'


# Convert results file from FastANI to a somewhat better format.
# This is the only rule whose output is needed in other files.
rule format_fastani_results:
	input: rules.fastani.output
	output: 'intermediate-data/fastani/{genomeset}-formatted.csv'
	params:
		filenames=lambda wc: get_genome_fasta_files(wc, full_path=False),
	script: '../scripts/format-fastani-pw-output.py'
