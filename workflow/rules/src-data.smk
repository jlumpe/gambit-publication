"""Download source data.

For rules that download multiple genomes one at a time, the files are placed in a separate directory
and then a symlink to that directory is created as the rule's output.
That way when some but not all downloads fail Snakemake will delete the symlink instead of the
directory containing the successful downloads.
The rule can then be rerun and it will not attempt to re-download the genomes already present.

These rules also have an nworkers parameter, which is the number of simultaneous downloads. This is
used instead of Snakemake's built-in "threads" variable because the download threads are IO-bound
instead of cpu-bound.
"""

# Store variations of resource files here for use in test mode
RESOURCES_TEST_DIR = 'resources/.test'

# Actually download genomes to this directory
GENOMES_DL_DIR = 'resources/genomes/.download'

# URL prefix for NCBI files
_ncbi_protocol = 'http' if config['src_data']['ncbi']['use_http'] else 'ftp'
NCBI_FTP_PREFIX = _ncbi_protocol + '://ftp.ncbi.nlm.nih.gov/'

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


def download_and_link(items, dl_path, link_path, nworkers, **kw):
	"""Download list of files to a directory and create a symlink to that directory."""

	from gambit_pub.download import download_parallel
	from gambit_pub.utils import symlink_to_relative

	dl_path = Path(dl_path)
	dl_path.mkdir(parents=True, exist_ok=True)
	download_parallel(items, dl_path, nworkers=nworkers, **kw)
	symlink_to_relative(dl_path, link_path)


def _get_genomeset(wc_or_gset):
	"""Get genome set ID string from argument."""
	if isinstance(wc_or_gset, str):  # Given the ID itself
		return wc_or_gset
	elif wc_or_gset is None:  # None, return wildcard placeholder
		return '{genomeset}'
	else:  # Assume wildcards object
		return wc_or_gset.genomeset

def get_genomes_list_file(wildcards_or_genomeset, test=TEST):
	"""Get the genomes.txt file for given genome set, using reduced version in test mode."""
	gset = _get_genomeset(wildcards_or_genomeset)
	parent_dir = RESOURCES_TEST_DIR if test else 'resources'
	return f'{parent_dir}/genomes/{gset}/genomes.txt'

def get_genomes_table_file(wildcards_or_genomeset, test=TEST):
	"""Get the genomes.txt file for given genome set, using reduced version in test mode."""
	gset = _get_genomeset(wildcards_or_genomeset)
	parent_dir = RESOURCES_TEST_DIR if test else 'resources'
	return f'{parent_dir}/genomes/{gset}/genomes.csv'


# Create truncated versions of genomes.txt when in test mode
rule truncated_genome_list:
	input: get_genomes_list_file(None, test=False)
	output: get_genomes_list_file(None, test=True)
	params:
		nlines=config['test_mode']['genome_cap'],
	shell:
		"head -n {params[n]} {input} > {output}"


# Create truncated versions of genomes.csv when in test mode
rule truncated_genome_table:
	input: get_genomes_table_file(None, test=False)
	output: get_genomes_table_file(None, test=True)
	params:
	      nlines=config['test_mode']['genome_cap'] + 1,
	shell:
		"head -n {params[n]} {input} > {output}"


# Download FASTA files for genome set 1 or 2 (both from NCBI FTP server)
rule get_genome_set_12:
	input: get_genomes_table_file
	output: directory('resources/genomes/{genomeset}/fasta/')
	params:
		dl_dir=f'{GENOMES_DL_DIR}/{{genomeset}}/fasta/',
		nworkers=config['src_data']['nworkers'],
		show_progress=config['show_progress'],  # Show progress bar
	wildcard_constraints:
		genomeset="set[12]",
	run:
		table = pd.read_csv(input[0])
		items = [
			(NCBI_FTP_PREFIX + row.ftp_path, row.assembly_accession + '.fa.gz', row.md5)
			for _, row in table.iterrows()
		]
		download_and_link(items, params['dl_dir'], output[0], params['nworkers'],
		                  progress=params['show_progress'], desc=f'{wildcards.genomeset} FASTA')


# Download FASTA files for genome sets 3 and 4
rule get_genome_set_34:
	input:
		get_genomes_list_file
	output:
		directory('resources/genomes/{genomeset}/fasta/')
	params:
		dl_dir=GENOMES_DL_DIR + '{genomeset}/fasta/',
		nworkers=config['src_data']['nworkers'],
		show_progress=config['show_progress'],  # Show progress bar
	wildcard_constraints:
		genomeset="set[34]",
	run:
		from gambit.util.io import read_lines

		gset = wildcards.genomeset
		gs_dir = config['src_data']['genome_sets'][gset]['fasta'].rstrip('/')
		prefix = GCS_PREFIX + gs_dir + '/'
		items = [(prefix + fname, fname, None) for fname in read_lines(input[0])]

		download_and_link(items, params['dl_dir'], output[0], params['nworkers'],
		                  progress=params['show_progress'], desc=f'{gset} FASTA')


# Download FASTQ files for genome set 3
rule get_genome_set_3_fastq:
	input:
		get_genomes_list_file('set3')
	output:
		directory('resources/genomes/set3/fastq/')
	params:
		dl_dir=GENOMES_DL_DIR + 'set3/fastq/',
		nworkers=config['src_data']['nworkers'],
		show_progress=config['show_progress'],  # Show progress bar
	run:
		gs_dir = config['src_data']['genome_sets']['set3']['fastq'].rstrip('/')
		prefix = GCS_PREFIX + gs_dir + '/'

		items = []
		with open(input[0]) as f:
			for line in f:
				fname = line.strip().rsplit('.', 1)[0] + '.fastq.gz'
				items.append((prefix + fname, fname, None))

		download_and_link(items, params['dl_dir'], output[0], params['nworkers'],
		                  progress=params['show_progress'], desc='set3 FASTQ')


# Download FASTA files for genome set 5
rule get_genome_set_5:
	output:
		directory('resources/genomes/set5/fasta'),
	params:
		url=GCS_PREFIX + config['src_data']['genome_sets']['set5']['tarball'],
	shell:
		"""
		mkdir {output}
		curl {params[url]} | tar -xzf - -C {output} --strip-components=1
		"""


# Download all source data
rule get_src_data:
	input:
		*rules.get_gambit_db.output,
		*expand('resources/genomes/{gset}/fasta', gset=['set1', 'set2', 'set3', 'set4', 'set5']),
		*rules.get_genome_set_3_fastq.output,
