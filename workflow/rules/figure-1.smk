rule figure_1:
	params:
		genome_sets=COMPARISON_GENOME_SETS,
	input:
		expand(rules.gambit_vs_ani.output, genomeset=COMPARISON_GENOME_SETS, k=K, prefix=PREFIX),
	output:
		figure='results/figure-1/figure-1.png',
		stats='results/figure-1/stats.csv',
	script: '../scripts/figure-1.py'
