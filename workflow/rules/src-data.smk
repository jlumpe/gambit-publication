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
rule fetch_gambit_db:
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


def _get_genomeset(wc_or_gset, require=False):
	"""Get genome set ID string from argument."""
	if isinstance(wc_or_gset, str):  # Given the ID itself
		return wc_or_gset
	elif wc_or_gset is None:  # None, return wildcard placeholder
		if require:
			raise TypeError('A value is required')
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

def get_genomes_fasta_files(wildcards_or_genomeset, test=TEST, full_path=True):
	"""Get paths of FASTA files for the given genome set (relative to root directory)."""
	gset = _get_genomeset(wildcards_or_genomeset, require=True)
	parent_dir = f'resources/genomes/{gset}/fasta/' if full_path else ''
	list_file = get_genomes_list_file(gset, test=test)
	return [os.path.join(parent_dir, filename) for filename in read_lines(list_file)]


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


def _get_gset_12_dl_items(gset):
	table = pd.read_csv(get_genomes_table_file(gset))
	return [
		(NCBI_FTP_PREFIX + row.ftp_path, row.assembly_accession + '.fa.gz', row.md5)
		for _, row in table.iterrows()
	]

def _get_gset_34_dl_items(gset):
	from gambit.util.io import read_lines

	gs_dir = config['src_data']['genome_sets'][gset]['fasta'].rstrip('/')
	prefix = GCS_PREFIX + gs_dir + '/'
	list_file = get_genomes_list_file(gset)
	return [(prefix + fname, fname, None) for fname in read_lines(list_file)]

# Download FASTA files for genome sets 1-4
rule fetch_genome_set_1234:
	output: directory('resources/genomes/{genomeset}/fasta')
	params:
		dl_dir=f'{GENOMES_DL_DIR}/{{genomeset}}/fasta/',
		nworkers=config['src_data']['nworkers'],
		show_progress=config['show_progress'],
	wildcard_constraints:
		genomeset="set[1234]",
	run:
		from gambit_pub.download import download_parallel
		from gambit_pub.utils import symlink_to_relative

		gset = wildcards.genomeset

		if gset in ('set1', 'set2'):
			items = _get_gset_12_dl_items(gset)
		elif gset in ('set3', 'set4'):
			items = _get_gset_34_dl_items(gset)
		else:
			raise ValueError(gset)

		dl_path = Path(params['dl_dir'])
		dl_path.mkdir(parents=True, exist_ok=True)
		download_parallel(items, dl_path, nworkers=params['nworkers'], progress=params['show_progress'])
		symlink_to_relative(dl_path, output[0])


# Download FASTQ files for genome set 3
# Unlike the FASTA download rules, this one operates on one FASTQ file at a time because they are
# large and not every one is used.
rule fetch_genome_set_3_fastq:
	output:
		'resources/genomes/set3/fastq/{genome}.fastq.gz'
	params:
		gs_dir=GCS_PREFIX + config['src_data']['genome_sets']['set3']['fastq'].rstrip('/'),
	shell:
		"curl {params[gs_dir]}/{wildcards[genome]}.fastq.gz -o {output}"


def get_genome_set_3_fastq_files(wildcards=None):
	list_file = get_genomes_list_file('set3')
	genomes = [file.rsplit('.', 1)[0] for file in read_lines(list_file)]
	return expand(rules.fetch_genome_set_3_fastq.output, genome=genomes)

# Download FASTQ files for all set 3 genomes
# Not needed by main figures and tables and not included in fetch_src_data rule
rule fetch_genome_set_3_fastq_all:
	input: get_genome_set_3_fastq_files


# Download FASTA files for genome set 5
rule fetch_genome_set_5:
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
rule fetch_src_data:
	input:
		*rules.fetch_gambit_db.output,
		# *expand(rules.fetch_genome_set_1234.output, genomeset=['set1', 'set2', 'set3', 'set4']),
		*expand(rules.fetch_genome_set_1234.output, genomeset=['set1', 'set2', 'set3']),
		*rules.fetch_genome_set_5.output,
		# TODO - only the needed Set 3 FASTQ files
