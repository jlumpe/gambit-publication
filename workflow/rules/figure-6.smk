rule figure_6:
	input:
	     genomes_csv=get_genomes_table_file('figure_6'),
	     pw_dists=expand(rules.gambit_pw_dists.output[0], genomeset='figure_6', k=K, prefix=PREFIX)[0],
	output:
	      "results/figure-6/figure-6.png"
	script:
	      "../scripts/figure-6.py"
