# Main 3 genome sets that are used for figures 1 and 2
COMPARISON_GENOME_SETS = ['ondov_2016', 'konstantinidis_2005']


### Figure 1 ###

rule figure_1:
	input:
		expand(rules.gambit_vs_ani.output[0], genomeset=COMPARISON_GENOME_SETS, k=K, prefix=PREFIX),
	# output:
	# 	"results/figure-1/figure-1.png"


### Figure 6 ###

rule figure_6:
	input:
		genomes_csv="resources/genomes/figure_6/genomes.csv",
		pw_dists=expand(rules.gambit_pw_dists.output[0], genomeset='figure_6', k=K, prefix=PREFIX)[0],
	output:
		"results/figure-6/figure-6.png"
	script:
	    "../scripts/figure-6.py"
