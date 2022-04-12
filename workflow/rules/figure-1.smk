rule figure_1:
	input:
		expand(rules.gambit_vs_ani.output[0], genomeset=COMPARISON_GENOME_SETS, k=K, prefix=PREFIX),
	output:
		touch('results/figure-1/done')  # TODO
