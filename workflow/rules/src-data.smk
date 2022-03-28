"""Download source data.

For rules that download multiple genomes one at a time, the files are placed in a separate directory
and then a symlink to that directory is created as the rule's output.
That way when some but not all downloads fail Snakemake will delete the symlink instead of the
directory containing the successful downloads.
The rule can then be rerun and it will not attempt to re-download the genomes already present.
"""

GENOMES_DL_DIR = 'resources/genomes/.download/'

# Prefix for GCS URLs
GCS_PREFIX = 'https://storage.googleapis.com/'


# Download GAMBIT database files
rule get_gambit_db:
	output:
		genomes='resources/gambit-db/db-genomes.db',
		signatures='resources/gambit-db/db-signatures.h5',
	params:
		genomes_url=GCS_PREFIX + config['src_data']['gambit_db']['genomes'],
		signatures_url=GCS_PREFIX + config['src_data']['gambit_db']['signatures'],
	shell:
	    '''
		curl -f -o {output.genomes} {params[genomes_url]}
		curl -f -o {output.signatures} {params[signatures_url]}
		'''


def download_and_link(items, dl_path, link_path, threads):
	"""Download list of files to a directory and create a symlink to that directory."""

	from gambit_pub.download import download_parallel
	from gambit_pub.utils import symlink_to_relative

	# Try for a multiple of what Snakemake tells us to use because the task is io/network-bound and
	# not cpu-bound, but enforce a maximum.
	nworkers = min(threads * 2, 12)

	dl_path = Path(dl_path)
	dl_path.mkdir(parents=True, exist_ok=True)
	download_parallel(items, dl_path, nworkers=nworkers)
	symlink_to_relative(dl_path, link_path)


# Creates truncated versions of genomes.csv and genomes.txt when in test mode
rule truncated_genome_list:
	output:
		"resources/genomes/{genomeset}/genomes-test.{ext}"
	wildcard_constraints:
		ext='txt|csv',
	run:
		out_dir = Path(output[0]).parent
		src = out_dir / ('genomes.' + wildcards['ext'])
		n = config['genome_cap']
		if wildcards['ext'] == 'csv':
			n += 1  # Account for header

		with open(src) as fsrc, open(output[0], 'w') as fdst:
			for i, line in enumerate(fsrc):
				if i >= n:
					break
				fdst.write(line)


# Get file containing list of genomes for the given genome set
def get_genomes_list_file(genomeset):
	fname = 'genomes-test.txt' if TEST else 'genomes.txt'
	return f'resources/genomes/{genomeset}/{fname}'

def genomes_list_file(wildcards):
	return get_genomes_list_file(wildcards.genomeset)

# Get file containing table of genome attributes for the given genome set
def get_genomes_table_file(genomeset):
	fname = 'genomes-test.csv' if TEST else 'genomes.csv'
	return f'resources/genomes/{genomeset}/{fname}'

def genomes_table_file(wildcards):
	return get_genomes_table_file(wildcards.genomeset)


# Download genome set 1 or 2 (both from NCBI FTP server)
rule get_genome_set_12:
	input:
		genomes_table_file
	output:
		directory('resources/genomes/{genomeset}/fasta/')
	params:
		dl_dir=GENOMES_DL_DIR + '{genomeset}/fasta/',
	wildcard_constraints:
		genomeset="set[12]",
	run:
		table = pd.read_csv(input[0])
		items = [(row.url, row.assembly_accession + '.fa.gz', row.md5) for _, row in table.iterrows()]
		download_and_link(items, params['dl_dir'], output[0], threads)


# Download genome set 3
rule get_genome_set_3:
	input:
		get_genomes_list_file('set3')
	output:
		directory('resources/genomes/set3/fasta/')
	params:
		dl_dir=GENOMES_DL_DIR + 'set3/fasta/',
	run:
		gs_dir = config['src_data']['genome_sets']['set3']['fasta'].rstrip('/')
		prefix = GCS_PREFIX + gs_dir + '/'

		items = []
		with open(input[0]) as f:
			for line in f:
				fname = line.strip()
				items.append((prefix + fname, fname, None))

		download_and_link(items, params['dl_dir'], output[0], threads)


# Download fastq files for genome set 3
rule get_genome_set_3_fastq:
	input:
		get_genomes_list_file('set3')
	output:
		directory('resources/genomes/set3/fastq/')
	params:
		dl_dir=GENOMES_DL_DIR + 'set3/fastq/',
	run:
		gs_dir = config['src_data']['genome_sets']['set3']['fastq'].rstrip('/')
		prefix = GCS_PREFIX + gs_dir + '/'

		items = []
		with open(input[0]) as f:
			for line in f:
				fname = line.strip().rsplit('.', 1)[0] + '.fastq.gz'
				items.append((prefix + fname, fname, None))

		download_and_link(items, params['dl_dir'], output[0], threads)


# Download genomes for figure 6
rule get_genomes_fig6:
	output:
	    directory('resources/genomes/figure_6/fasta'),
	params:
		url=GCS_PREFIX + config['src_data']['genome_sets']['figure_6']['tarball'],
	shell:
		"""
		parent=$(dirname {output})
		curl {params[url]} | tar -xzf - -C $parent
		"""


# Download all source data
rule get_src_data:
	input:
		# *rules.get_gambit_db.output,
		*expand('resources/genomes/{gset}/fasta', gset=['set1', 'set2', 'set3', 'figure_6']),
		# rules.get_genome_set_3_fastq.output,
