"""
Generate all figures and some related data tables.
"""


### Figure 1 ###

rule figure_1:
	params:
		genome_sets=COMPARISON_GENOME_SETS,
	input:
		expand(rules.gambit_vs_ani.output, genomeset=COMPARISON_GENOME_SETS, k=K, prefix=PREFIX),
	output:
		figure='results/figure-1/figure-1.png',
		stats='results/figure-1/stats.csv',
	script: '../scripts/figure-1.py'


### Figure 2 and supplemental figure 1 ###

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
# Paramspace wildcard must be one of several predefined strings.
rule gambit_ani_correlation:
	params:
	      params_df=gambit_paramspace,
	input: gambit_ani_correlation_input,
	output: 'results/gambit-ani-correlation/{paramspace}.csv'
	script: '../scripts/gambit-vs-ani-correlation.py'


rule figure_2a:
	input: expand(rules.gambit_ani_correlation.output, paramspace='prefix_length')
	output: 'results/figure-2/figure-2a.png'
	run:
		import gambit_pub.paramspace_exploration as pex

		paramdata = pex.get_param_data(input[0], config)

		pex.set_style()
		fg = pex.spearman_vs_k(
			paramdata,
			col='prefix_len',
		)

		for plen, ax in fg.axes_dict.items():
			ax.set_title(paramdata.prefix_map[plen, 0])

		pex.highlight_default_axis(fg.axes_dict[paramdata.dflt_plen])

		fg.figure.savefig(output[0])

rule figure_2b:
	input: expand(rules.gambit_ani_correlation.output, paramspace='prefix_sequence')
	output: 'results/figure-2/figure-2b.png'
	run:
		import gambit_pub.paramspace_exploration as pex

		paramdata = pex.get_param_data(input[0], config)

		pex.set_style()
		fg = pex.spearman_vs_k(
			paramdata,
			col='prefix_version',
			col_wrap=4,
		)

		for pver, ax in fg.axes_dict.items():
			ax.set_title(paramdata.prefix_map[paramdata.dflt_plen, pver])

		pex.highlight_default_axis(fg.axes_dict[paramdata.dflt_pver])

		fg.figure.savefig(output[0])

rule figure_2:
	input:
		rules.figure_2a.output,
		rules.figure_2b.output,
	output: touch('results/figure-2/.completed')


rule supplemental_figure_1:
	input: expand(rules.gambit_ani_correlation.output, paramspace='full')
	output: 'results/supplemental-figure-1/supplemental-figure-1.png'
	run:
		import gambit_pub.paramspace_exploration as pex
		import matplotlib as mpl

		paramdata = pex.get_param_data(input[0], config)

		pex.set_style()
		mpl.rcParams.update({
			'axes.titlesize': 12,
		})

		fg = pex.spearman_vs_k(
			paramdata,
			col='prefix_version',
			row='prefix_len',
			height=2,
		)

		for (plen, pver), ax in fg.axes_dict.items():
			ax.set_title(paramdata.prefix_map[plen, pver])

		pex.highlight_default_axis(fg.axes_dict[paramdata.dflt_plen, paramdata.dflt_pver])

		fg.figure.savefig(output[0])


### Figure 3 ###

# Compare signatures derived from single FASTQ file to standard FASTA
rule fastq_signatures:
	input:
		signatures=expand(rules.gambit_signatures.output, genomeset='set3', k=K, prefix=PREFIX)[0],
		fastq='resources/genomes/set3/fastq/{file}.fastq.gz',
	output:
		'intermediate-data/fastq-signatures/{file}.csv'
	params:
		fasta_name='{file}.fasta',
	wildcard_constraints:
		file=r'[\w-]+',
	script:
		'../scripts/fastq-signatures.py'


def get_fig3_input(wildcards):
	with open(get_genomes_list_file('set3')) as f:
		names = [line.split('.')[0] for line in f]

	return expand(rules.fastq_signatures.output, file=names)


rule figure_3:
	input:
		get_fig3_input
	output:
		touch('results/figure-3/done')  # TODO


### Figure 6 ###

rule figure_6:
	input:
		genomes_csv=get_genomes_table_file('set5'),
		pw_dists=expand(rules.gambit_pw_dists.output[0], genomeset='set5', k=K, prefix=PREFIX)[0],
	output:
		"results/figure-6/figure-6.png"
	script:
		"../scripts/figure-6.py"
