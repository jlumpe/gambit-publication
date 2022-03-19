

# Generate signatures for genomes in each genome set
rule genome_set_signatures:
	input:
		"resources/genomes/{genomeset}/fasta",
	output:
		"intermediate-data/genomes/{genomeset}/signatures/{k}-{prefix}.h5",
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


# Calculate pairwise distances from genome set signatures
rule genome_set_pw_dists:
	input:
		lambda wildcards: expand(rules.genome_set_signatures.output[0], **wildcards),
	output:
		"intermediate-data/genomes/{genomeset}/pw-dists/{k}-{prefix}.csv",
	shell:
	    "gambit dist -s --qs {input[0]} -o {output[0]}"
