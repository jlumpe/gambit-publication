"""
Benchmark performance of GAMBIT and other tools.
"""



def get_dist_benchmark_input(wildcards):
	genome_params = config['dist_benchmarks']['genomes'][wildcards.genomes_key]
	query_set = genome_params['queries']['genome_set']
	ref_set = genome_params['refs']['genome_set']
	return dict(
		queries_fasta_dir=get_genomes_fasta_dir(query_set),
		queries_list_file=get_genomes_list_file(query_set),
		refs_fasta_dir=get_genomes_fasta_dir(ref_set),
		refs_list_file=get_genomes_list_file(ref_set),
	)

# Run distance calculation benchmarks for a particular set of query/reference genomes defined in config
rule benchmark_dists:
	input: unpack(get_dist_benchmark_input)
	output:
		table='results/benchmarks/dists/{genomes_key}.csv',
		extra='results/benchmarks/dists/{genomes_key}.json',
	wildcard_constraints:
		genomes=r'\w+',
	# This is just to prevent other jobs from running at the same time, doesn't influence the number
	# of threads/cores used in benchmarks.
	threads: workflow.cores
	shadow: 'minimal'
	script: '../scripts/benchmark-dists.py'


# Run distance benchmarks for all sets of query/reference genomes defined in config
rule benchmark_dists_all:
	input: expand(rules.benchmark_dists.output, genomes_key=list(config['dist_benchmarks']['genomes']))
	output: touch('results/benchmarks/dists/.completed')
