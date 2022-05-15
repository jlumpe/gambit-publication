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
		genomes_csv=get_genomes_table_file('set4'),
		pw_dists=expand(rules.gambit_pw_dists.output[0], genomeset='set4', k=K, prefix=PREFIX)[0],
	output:
		"results/figure-6/figure-6.png"
	script:
		"../scripts/figure-6.py"
