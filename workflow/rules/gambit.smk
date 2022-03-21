"""Run GAMBIT CLI tool."""


# Generate signatures for genomes in each genome set
rule gambit_signatures:
	input:
		"resources/genomes/{genomeset}/fasta",
	output:
		"intermediate-data/signatures/{genomeset}-{k}-{prefix}.h5",
	shell:
		"""
		genomes_dir=$(dirname {input[0]})
		gambit signatures create \
			-k {wildcards.k} \
			-p {wildcards.prefix} \
			-l $genomes_dir/genomes.txt \
			--ldir {input[0]} \
			-o {output[0]}
		"""


# Pairwise GAMBIT distances from genome set signatures
rule gambit_pw_dists:
	input:
		rules.gambit_signatures.output
	output:
		"intermediate-data/gambit-pw-dists/{genomeset}-{k}-{prefix}.csv",
	shell:
	    "gambit dist -s --qs {input[0]} -o {output[0]}"
