"""
Generate all figures and some related data tables.
"""


### Figure 1 ###

rule figure_1:
	params:
		genome_sets=COMPARISON_GENOME_SETS,
	input:
		pairs=expand(rules.gambit_vs_ani.output['pairs'], genomeset=COMPARISON_GENOME_SETS, k=K, prefix=PREFIX),
		stats=expand(rules.gambit_vs_ani.output['stats'], genomeset=COMPARISON_GENOME_SETS, k=K, prefix=PREFIX),
	output:
		figure='results/figures/figure-1.png',
		stats='results/figures/figure-1.csv',
	script: '../scripts/figure-1.py'


### Figure 2 ###

rule figure_2a:
	input: expand(rules.gambit_ani_correlation.output, paramspace='prefix_length')
	output: 'results/figures/figure-2a.png'
	params:
		paramspace='prefix_length',
	script: '../scripts/figure-2.py'

rule figure_2b:
	input: expand(rules.gambit_ani_correlation.output, paramspace='prefix_sequence')
	output: 'results/figures/figure-2b.png'
	params:
		paramspace='prefix_sequence',
	script: '../scripts/figure-2.py'

rule figure_2:
	input:
		rules.figure_2a.output,
		rules.figure_2b.output,
	output: touch('results/figures/.figure-2-completed')


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
	output: 'results/figures/figure-3.png'
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
	output: 'results/figures/figure-4{subplot}.png'
	wildcard_constraints:
		subplot='[a-z]'
	threads: workflow.cores
	params:
		conf=lambda wildcards: config['figure_4']['subplots'][wildcards.subplot],
	script: '../scripts/figure-4-subplot.py'

# Generate all subplots of figure 4
rule figure_4:
	input: expand(rules.figure_4_subplot.output, subplot=list(config['figure_4']['subplots']))
	output: touch('results/figures/.figure-4-completed')


### Figure 5 ###

# Generate one of figure 5's subplots
# List of subplot parameters is in config file under figure_5/subplots
rule figure_5_subplot:
	input:
		db_signatures=rules.fetch_gambit_db.output['signatures'],
		db_genomes=rules.fetch_gambit_db.output['genomes'],
	output: 'results/figures/figure-5{subplot}.png'
	wildcard_constraints:
		subplot='[a-z]'
	threads: workflow.cores
	params:
		conf=lambda wildcards: config['figure_5']['subplots'][wildcards.subplot],
	script: '../scripts/figure-5-subplot.py'

# Generate all subplots of figure 5
rule figure_5:
	input: expand(rules.figure_5_subplot.output, subplot=list(config['figure_5']['subplots']))
	output: touch('results/figures/.figure-5-completed')


### Figure 6 ###

rule figure_6:
	input:
		genomes_csv=get_genomes_table_file('set5'),
		pw_dists=expand(rules.gambit_pw_dists.output[0], genomeset='set5', k=K, prefix=PREFIX)[0],
	output:
		"results/figures/figure-6.png"
	script:
		"../scripts/figure-6.py"
