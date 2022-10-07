"""Run GAMBIT CLI tool."""


# Generate signatures for genomes in each genome set
rule gambit_signatures:
	input:
		fasta_dir=get_genomes_fasta_dir,
		list_file=get_genomes_list_file,
	output:
		"intermediate-data/signatures/{genomeset}-{k}-{prefix}.gs",
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


# Comparison of GAMBIT distance and ANI (via FastANI)
rule gambit_vs_ani:
	input:
		gambit=rules.gambit_pw_dists.output[0],
		fastani=rules.format_fastani_results.output[0],
	output:
		pairs='intermediate-data/gambit-vs-ani/{genomeset}-{k}-{prefix}.csv',
		stats='intermediate-data/gambit-vs-ani/{genomeset}-{k}-{prefix}.json',
	script: '../scripts/gambit-vs-ani.py'


@cache
def get_gambit_paramspace(wildcards):
	"""Parameter space to explore for GAMBIT distance / ANI correlation.

	The "paramspace" wildcard must be one of several predefined options.
	"""
	values = config['gambit_param_space']

	if wildcards.paramspace == 'full':
		k_vals = values['k']
		prefix_lens = values['prefix_len']
		base_prefixes = values['base_prefix']
	elif wildcards.paramspace == 'prefix_length':
		k_vals = values['k']
		prefix_lens = values['prefix_len']
		base_prefixes = [values['base_prefix'][0]]
	elif wildcards.paramspace == 'prefix_sequence':
		k_vals = values['k']
		prefix_lens = [len(PREFIX)]
		base_prefixes = values['base_prefix']
	elif wildcards.paramspace == 'default_only':
		k_vals = [K]
		prefix_lens = [len(PREFIX)]
		base_prefixes = [values['base_prefix'][0]]
	else:
		raise ValueError(f'Invalid value for "paramspace" wildcard: {wildcards.paramspace}')

	index = pd.MultiIndex.from_product(
		[COMPARISON_GENOME_SETS, k_vals, prefix_lens, range(len(base_prefixes))],
		names=['genomeset', 'k', 'prefix_len', 'prefix_version'],
	)
	df = index.to_frame(False)

	df['base_prefix'] = [base_prefixes[i] for i in df['prefix_version']]
	df['prefix'] = [row.base_prefix[:row.prefix_len] for _, row in df.iterrows()]

	return df


def gambit_ani_correlation_input(wildcards):
	"""Input for gambit_ani_correlation rule."""
	params_df = get_gambit_paramspace(wildcards)
	return [
		expand(rules.gambit_vs_ani.output['stats'], **row.to_dict())[0]
		for _, row in params_df.iterrows()
	]

# Consolidate GAMBIT distance vs ANI statistics for a range of parameter values.
# "paramspace" wildcard must be one of several predefined strings.
rule gambit_ani_correlation:
	input: gambit_ani_correlation_input,
	output: 'results/gambit-ani-correlation/{paramspace}.csv'
	params:
		params_df=get_gambit_paramspace,
	run:
		import json, pandas as pd

		stats_rows = []
		for stats_file in input:
			with open(stats_file) as f:
				stats_rows.append(json.load(f))

		stats_df = pd.DataFrame(stats_rows)

		combined = pd.concat([params['params_df'], stats_df], axis=1)
		combined.to_csv(output[0], index=False)
