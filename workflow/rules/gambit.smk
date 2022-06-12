"""Run GAMBIT CLI tool."""


# Generate signatures for genomes in each genome set
rule gambit_signatures:
	input:
		fasta_dir=get_genomes_fasta_dir,
		list_file=get_genomes_list_file,
	output:
		"intermediate-data/signatures/{genomeset}-{k}-{prefix}.h5",
	threads: workflow.cores
	shell:
		"""
		gambit signatures create \
			-c {threads} \
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
		"""
		gambit dist -s \
			-c {threads} \
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


@cache
def gambit_paramspace(wildcards):
	"""Parameter space to explore for GAMBIT distance / ANI correlation.

	The "paramspace" wildcard must be one of several predefined options.
	"""
	from gambit_pub.paramspace_exploration import make_params_df

	values = config['gambit_param_space']

	if wildcards.paramspace == 'full':
		args = (values['k'], values['prefix_len'], values['base_prefix'])
	elif wildcards.paramspace == 'prefix_length':
		args = (values['k'], values['prefix_len'], values['base_prefix'][0])
	elif wildcards.paramspace == 'prefix_sequence':
		args = (values['k'], len(PREFIX), values['base_prefix'])
	# elif wildcards.paramspace == 'default_only':
	# 	args = (K, len(PREFIX), PREFIX)
	else:
		raise ValueError(f'Invalid value for "paramspace" wildcard: {wildcards.paramspace}')

	return make_params_df(*args, COMPARISON_GENOME_SETS)

def gambit_ani_correlation_input(wildcards):
	"""Input for gambit_ani_correlation rule."""
	params_df = gambit_paramspace(wildcards)
	return [
		expand(rules.gambit_vs_ani.output, **row.to_dict())[0]
		for _, row in params_df.iterrows()
	]

# Calculate spearman correlation of GAMBIT distance vs ANI for a range of parameter values.
# "paramspace" wildcard must be one of several predefined strings.
rule gambit_ani_correlation:
	params:
		params_df=gambit_paramspace,
	input: gambit_ani_correlation_input,
	output: 'results/gambit-ani-correlation/{paramspace}.csv'
	script: '../scripts/gambit-vs-ani-correlation.py'
