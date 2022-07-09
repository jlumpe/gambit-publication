"""
Generate supplemental figures and data.
"""


rule sfig1:
	input: expand(rules.gambit_ani_correlation.output, paramspace='full')
	output: 'results/figures/supplemental-figure-1.png'
	params:
		paramspace='full',
	script: '../scripts/figure-2.py'


rule sfig2:
	input: expand(rules.gambit_vs_ani.output['pairs'], k=K, prefix=PREFIX, genomeset=COMPARISON_GENOME_SETS)
	output:
		figure='results/figures/supplemental-figure-2.png',
	params:
		genome_sets=COMPARISON_GENOME_SETS,
	script: '../scripts/supplemental-figure-2.py'
