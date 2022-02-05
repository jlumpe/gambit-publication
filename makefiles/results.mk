# Rules for running notebooks and generating results

.PHONY: results


# Generate all results
results: results/figure6


results/figure-6:
	# Generate figure 6
	$(conda_run) $(run_notebook) notebooks/figure-6/figure-6.ipynb
