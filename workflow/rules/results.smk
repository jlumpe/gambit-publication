# Generate figure 6
rule figure_6:
	input:
		genomes_csv="src-data/genomes/figure-6/genomes.csv",
		pw_dists=expand(rules.genome_set_pw_dists.output[0], genomeset='figure-6', k=K, prefix=PREFIX)[0],
	output:
		"results/figure-6.png"
	script:
	    "../scripts/figure-6.py"
