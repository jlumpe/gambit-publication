# Rules for running notebooks and generating results

.PHONY: results


# Generate all results
results: results/figure6


# Figure 6
results/figure-6:
	$(CONDA_RUN) $(RUN_NOTEBOOK) notebooks/figure-6/figure-6.ipynb
