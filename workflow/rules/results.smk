

# Generate figure 6
rule figure_6:
	input:
		rules.get_genomes_fig6.output[0]
	output:
		directory("results/figure-6")
	notebook:
		"../notebooks/figure-6.ipynb"
