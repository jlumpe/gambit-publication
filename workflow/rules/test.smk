"""
Test mode stuff.
"""


# Create truncated versions of genomes.txt when in test mode
checkpoint truncated_genome_list:
	input: lambda wc: get_genomes_list_file(wc, test=False)
	output: get_genomes_list_file(None, test=True)
	params:
		n=config['test_genome_cap'],
	shell:
		"head -n {params[n]} {input} > {output}"

# Create truncated versions of genomes.csv when in test mode
checkpoint truncated_genome_table:
	input: lambda wc: get_genomes_table_file(wc, test=False)
	output: get_genomes_table_file(None, test=True)
	params:
		n=config['test_genome_cap'] + 1,
	shell:
		"head -n {params[n]} {input} > {output}"
