"""
Generate supplemental figures and data.
"""


### Supplemental Figure 1 ###

rule supplemental_figure_1:
	input: expand(rules.gambit_ani_correlation.output, paramspace='full')
	output: 'results/supplemental-figure-1/supplemental-figure-1.png'
	run:
		import gambit_pub.paramspace_exploration as pex
		import matplotlib as mpl

		paramdata = pex.get_param_data(input[0], config)

		pex.set_style()
		mpl.rcParams.update({
			'axes.titlesize': 12,
		})

		fg = pex.spearman_vs_k(
			paramdata,
			col='prefix_version',
			row='prefix_len',
			height=2,
		)

		for (plen, pver), ax in fg.axes_dict.items():
			ax.set_title(paramdata.prefix_map[plen, pver])

		pex.highlight_default_axis(
			fg.axes_dict[paramdata.dflt_plen, paramdata.dflt_pver],
			paramdata.dflt_k,
		)

		fg.figure.savefig(output[0])
