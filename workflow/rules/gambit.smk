"""Run GAMBIT CLI tool."""


# Generate signatures for genomes in each genome set
rule gambit_signatures:
	input:
		fasta_dir=get_genomes_fasta_dir,
		list_file=get_genomes_list_file,
	output:
		"intermediate-data/signatures/{genomeset}-{k}-{prefix}.h5",
	# threads: workflow.cores  # This task is mostly IO-bound so let Snakemake run several at a time
	shell:
		"""
		gambit signatures create \
			-k {wildcards.k} \
			-p {wildcards.prefix} \
			-l {input[list_file]} \
			--ldir {input[fasta_dir]} \
			-o {output[0]} \
			--no-progress
		"""


# Pairwise GAMBIT distances from genome set signatures
rule gambit_pw_dists:
	input:
		rules.gambit_signatures.output
	output:
		"intermediate-data/gambit-pw-dists/{genomeset}-{k}-{prefix}.csv",
	threads: workflow.cores
	shell:
		"""gambit dist -s \
			--qs {input[0]} \
			-o {output[0]} \
			--no-progress
		"""


rule gambit_vs_ani:
	input:
		gambit=rules.gambit_pw_dists.output[0],
		fastani=rules.fastani.output[0],
	output:
		"intermediate-data/gambit-vs-ani/{genomeset}-{k}-{prefix}.csv",
	script:
		"../scripts/gambit-vs-ani.py"
