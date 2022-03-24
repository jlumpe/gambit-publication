"""
Compute all pairwise FastANI comparisons for a given genome set.

Genome set 1 is much larger than the others, and so needs to be divided up into several chunks
to make things reasonable. The size of each chunk (number of query genomes) is determined by the
genome_sets.{genomeset}.fastani_chunk_size config variable, if present.
"""


FASTANI_CHUNKS_DIR = "intermediate-data/fastani/chunks/{genomeset}"


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
		list_file="resources/genomes/{genomeset}/genomes.txt",
	output:
	      directory(f"{FASTANI_CHUNKS_DIR}/queries")
	run:
		out_dir = Path(output[0])
		out_dir.mkdir()

		# Ready query list file
		with open(input['list_file']) as f:
			lines = list(f)

		# Get chunk size from config (if present), default to single chunk
		gset_conf = config['genome_sets'].get(wildcards.genomeset, dict())
		chunk_size = gset_conf.get('fastani_chunk_size', len(lines))

		# Write chunked list files
		for begin in range(0, len(lines), chunk_size):
			chunk = lines[begin:begin + chunk_size]
			end = begin + len(chunk)

			with open(out_dir / f'{begin+1:03d}-{end:03d}.txt', 'w') as f:
				f.writelines(chunk)


# Generate single "chunk" of FastANI results, using all genomes for references but a subset of
# genomes for queries.
rule fastani_chunk:
	input:
		fasta="resources/genomes/{genomeset}/fasta",
		query_chunks=rules.fastani_query_chunks.output[0],
	output:
		temporary(f"{FASTANI_CHUNKS_DIR}/{{chunk}}.tsv")
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
		 cd {input[fasta]}
		 fastANI -k {params.k} --fragLen {params.fraglen} -t {threads} \\
			 --ql "$queries" --rl ../genomes.txt -o "$output"
		 """


# Get complete result file. This is the only rule whose output is needed in other files.
rule fastani:
	input:
	     get_fastani_input
	output:
	      "intermediate-data/fastani/{genomeset}.tsv"
	shell:
	     # Concatenate all chunks into a single output file
	     "cat {input} > {output}"
