"""
Rules that are not directly needed to generate the main or supplemental figures/data.
"""

# Run QUAST on assemblies in a genome set
# (QUAST isn't installed into environment by default)
rule genome_set_quast:
	input: 'resources/genomes/{genomeset}/fasta'
	output: directory('extra/quast/{genomeset}')
	threads: workflow.cores
	params:
		files=get_genomes_fasta_files,
	shell:
		"""
		quast --fast --silent -o {output} -t {threads} {params[files]}
		"""


# Generate genomes.csv for genome sets 3 and 4
# These are added to version control so it's not necessary to run normally
rule set_34_genomes_csv:
	input: rules.genome_set_quast.output
	output: 'resources/genomes/{genomeset}/genomes.csv'
	wildcard_constraints:
		genomeset='set[34]',
	params:
		filenames=lambda wc: get_genomes_fasta_files(wc, full_path=False),
	run:
		results = pd.read_csv(os.path.join(input[0], 'transposed_report.tsv'), sep='\t')

		# QUAST replaces dashses with underscores in file names, replace with original file names
		# but check consistency
		real_ids = [filename.split('.')[0] for filename in params['filenames']]
		quast_ids = results['Assembly']

		assert all(rid.replace('-', '_') == qid for rid, qid in zip(real_ids, quast_ids))

		results = results[['# contigs (>= 0 bp)', 'Total length (>= 0 bp)', 'N50', 'L50']]
		results.columns = ['n_contigs', 'total_length', 'N50', 'L50']
		results.index = pd.Series(real_ids, name='id')

		results.to_csv(output[0])


def get_fastq_kmers_all_input(wildcards=None):
	from gambit_pub.utils import stripext
	genomes = list(map(stripext, get_genomes_fasta_files('set3', full_path=False)))
	return expand(rules.fastq_kmers.output, genome=genomes)

# Run the fastq_kmers rule for all genomes in set 3
rule fastq_kmers_all:
	input: get_fastq_kmers_all_input
