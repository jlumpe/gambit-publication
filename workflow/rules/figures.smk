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


### Figure 2 ###

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

		pex.highlight_default_axis(fg.axes_dict[paramdata.dflt_plen], paramdata.dflt_k)

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

		pex.highlight_default_axis(fg.axes_dict[paramdata.dflt_pver], paramdata.dflt_k)

		fg.figure.savefig(output[0])

rule figure_2:
	input:
		rules.figure_2a.output,
		rules.figure_2b.output,
	output: touch('results/figure-2/.completed')


### Figure 3 ###

# Find k-mers in set3 FASTQ reads and compare to signatures derived from assembly FASTA
rule fastq_kmers:
	input:
		signatures=expand(rules.gambit_signatures.output, genomeset='set3', k=K, prefix=PREFIX)[0],
		fastq=rules.fetch_genome_set_3_fastq.output[0],
		genomes_table=get_genomes_table_file('set3'),
	output: directory('intermediate-data/fastq-kmers/{genome}/')
	params:
		fasta_id='{genome}',
		min_phred=[0, 10, 20, 30],
		truncate_reads=100_000 if TEST else None,
	wildcard_constraints:
		genome=r'[\w-]+',
	script:
		'../scripts/fastq-kmers.py'


def get_figure_3_genomes(wildcards=None):
	from gambit.cli.common import strip_seq_file_ext

	if TEST:
		fasta_files = get_genome_fasta_files('set3', test=True, full_path=False)
		return [strip_seq_file_ext(f) for f in fasta_files[:2]]
	else:
		return config['figure_3']['genomes']

def get_figure_3_input(wildcards=None):
	return expand(rules.fastq_kmers.output, genome=get_figure_3_genomes())

rule figure_3:
	input: get_figure_3_input
	output: 'results/figure-3/figure-3.png'
	params:
		genomes=get_figure_3_genomes,
		min_phred=20,
	script: '../scripts/figure-3.py'


### Figure 4 ###

# Generate one of figure 4's subplots
# List of subplot parameters is in config file under figure_4/subplots
rule figure_4_subplot:
	input:
		db_signatures=rules.fetch_gambit_db.output['signatures'],
		db_genomes=rules.fetch_gambit_db.output['genomes'],
	output: 'results/figure-4/figure-4{subplot}.png'
	wildcard_constraints:
		subplot='[a-z]'
	threads: workflow.cores
	params:
		conf=lambda wildcards: config['figure_4']['subplots'][wildcards.subplot],
	script: '../scripts/figure-4-subplot.py'

# Generate all subplots of figure 4
rule figure_4:
	input: expand(rules.figure_4_subplot.output, subplot=list(config['figure_4']['subplots']))
	output: touch('results/figure-4/.completed')


### Figure 5 ###

# Generate one of figure 5's subplots
# List of subplot parameters is in config file under figure_5/subplots
rule figure_5_subplot:
	input:
		db_signatures=rules.fetch_gambit_db.output['signatures'],
		db_genomes=rules.fetch_gambit_db.output['genomes'],
	output: 'results/figure-5/figure-5{subplot}.png'
	wildcard_constraints:
		subplot='[a-z]'
	threads: workflow.cores
	params:
		conf=lambda wildcards: config['figure_5']['subplots'][wildcards.subplot],
	script: '../scripts/figure-5-subplot.py'

# Generate all subplots of figure 5
rule figure_5:
	input: expand(rules.figure_5_subplot.output, subplot=list(config['figure_5']['subplots']))
	output: touch('results/figure-5/.completed')


### Figure 6 ###

rule figure_6:
	input:
		genomes_csv=get_genomes_table_file('set5'),
		pw_dists=expand(rules.gambit_pw_dists.output[0], genomeset='set5', k=K, prefix=PREFIX)[0],
	output:
		"results/figure-6/figure-6.png"
	script:
		"../scripts/figure-6.py"
