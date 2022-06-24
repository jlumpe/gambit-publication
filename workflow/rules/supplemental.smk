"""
Generate supplemental figures and data.
"""


### Supplemental Figure 1 ###

rule supplemental_figure_1:
	input: expand(rules.gambit_ani_correlation.output, paramspace='full')
	output: 'results/supplemental-figure-1/supplemental-figure-1.png'
	params:
		paramspace='full',
	script: '../scripts/figure-2.py'
