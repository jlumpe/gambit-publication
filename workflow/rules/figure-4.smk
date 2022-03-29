# Compare signatures derived from single FASTQ file to standard FASTA
rule fastq_signatures:
	input:
		signatures=expand(rules.gambit_signatures.output, genomeset='set3', k=K, prefix=PREFIX)[0],
		fastq='resources/genomes/set3/fastq/{file}.fastq.gz',
	output:
		'intermediate-data/fastq-signatures/{file}.csv'
	params:
		fasta_name='{file}.fasta',
	wildcard_constraints:
		file=r'[\w-]+',
	threads: 1  # Not parallelized
	script:
		'../scripts/fastq-signatures.py'


def get_fig4_input(wildcards):
	with open(get_genomes_list_file('set3')) as f:
		names = [line.split('.')[0] for line in f]

	names = [n for n in names if os.path.isfile(expand(rules.fastq_signatures.input['fastq'], file=n)[0])]
	files = expand(rules.fastq_signatures.output, file=names)
	return files


rule figure_4:
	input:
	     get_fig4_input
