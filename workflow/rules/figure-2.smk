def get_fig2_params():
	"""Data frame of all combinations of parameter values for figure 2."""

	values = config['figure_2']
	index = pd.MultiIndex.from_product(
		[COMPARISON_GENOME_SETS, values['k'], values['prefix_len'], range(len(values['base_prefix']))],
		names=['genomeset', 'k', 'prefix_len', 'prefix_version'],
	)
	df = index.to_frame(False)

	df['base_prefix'] = [values['base_prefix'][v] for v in df['prefix_version']]
	df['prefix'] = [row.base_prefix[:row.prefix_len] for _, row in df.iterrows()]
	df['input_key'] = [
		'{genomeset}-{k}-{prefix_len}-{prefix_version}'.format(**row.to_dict())
		for _, row in df.iterrows()
	]

	return df

def get_fig2_input(wildcards=None):
	params = get_fig2_params()
	return {
		row.input_key: expand(rules.gambit_vs_ani.output, **row.to_dict())[0]
		for _, row in params.iterrows()
	}

# Generate table of gambit/fastani correlations
rule figure_2_stats:
	params:
	      params_df=get_fig2_params(),
	input:
	     unpack(get_fig2_input)
	output:
	      stats="results/figure-2/gambit-ani-correlation.csv",
	script:
	      "../scripts/figure-2-stats.py"

rule figure_2:
	input:
	     rules.figure_2_stats.output
	output:
	      "results/figure-2/figure-2.png"
	script:
	      "../scripts/figure-2.py"
