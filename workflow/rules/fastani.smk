"""
Compute all pairwise FastANI comparisons for a given genome set.

The Ondov-2016 set is much larger than the others, and so needs to be divided up into several chunks
to make things reasonable. The size of each chunk (number of query genomes) is determined by the
genome_sets.{genomeset}.fastani_chunk_size config variable, if present.
"""


def fastani_input(wildcards):
	"""Get input files for the fastani rule (output chunks to combine).

	This depends on the output of the fastani_chunks checkpoint.
	"""
	chunkdir = Path(checkpoints.fastani_chunks.get(genomeset=wildcards.genomeset).output[0])
	chunks = [f.stem for f in chunkdir.glob('*-*.txt')]
	return expand(rules.fastani_chunk.output[0], genomeset=wildcards.genomeset, chunk=chunks)


# Divide query genome list into chunks so that FastANI run time for each is reasonable.
checkpoint fastani_chunks:
	output:
	      directory("intermediate-data/genomes/{genomeset}/fastani/chunks")
	run:
		out_dir = Path(output[0])
		out_dir.mkdir()

		# Ready query list file
		with open(f'resources/genomes/{wildcards.genomeset}/genomes.txt') as f:
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
		"resources/genomes/{genomeset}/fasta"
	output:
		temporary("intermediate-data/genomes/{genomeset}/fastani/fastani-{chunk}.tsv")
	wildcard_constraints:
		chunk="\d+-\d+"
	params:
		k=config['fastani']['k'],
		fraglen=config['fastani']['fraglen'],
	threads: workflow.cores
	shell:
	     """
		 output="$(realpath {output})"
		 queries="$(realpath intermediate-data/genomes/{wildcards.genomeset}/fastani/chunks/{wildcards.chunk}.txt)"
		 cd {input}
		 fastANI -k {params.k} --fragLen {params.fraglen} -t {threads} \\
			 --ql "$queries" --rl ../genomes.txt -o "$output"
		 """


# Get complete result file. This is the only rule whose output is needed in other files.
rule fastani:
	input:
	     fastani_input
	output:
	      "intermediate-data/genomes/{genomeset}/fastani/fastani.tsv"
	shell:
	     # Concatenate all chunks into a single output file
	     "cat {input} > {output}"
