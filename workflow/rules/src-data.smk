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

# Actually download genomes to this directory
GENOMES_DL_DIR = f'{DL_RESOURCES}/genomes/.download'

# URL prefix for NCBI files
_ncbi_protocol = 'http' if config['ncbi_use_http'] else 'ftp'
NCBI_FTP_PREFIX = _ncbi_protocol + '://ftp.ncbi.nlm.nih.gov/'

# Prefix for GCS URLs
GCS_PREFIX = 'https://storage.googleapis.com/'


# Download GAMBIT database files
rule fetch_gambit_db:
	output:
		genomes=f'{DL_RESOURCES}/gambit-db/db-genomes.db',
		signatures=f'{DL_RESOURCES}/gambit-db/db-signatures.h5',
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
	"""Get the genomes.txt file for given genome set, using truncated version in test mode."""
	gset = _get_genomeset(wildcards_or_genomeset)
	parent_dir = 'intermediate-data/test/genomes' if test else f'{SRC_DIR}/resources/genomes'
	return f'{parent_dir}/{gset}/genomes.txt'

def get_genomes_table_file(wildcards_or_genomeset, test=TEST):
	"""Get the genomes.csv file for given genome set, using truncated version in test mode."""
	gset = _get_genomeset(wildcards_or_genomeset)
	parent_dir = 'intermediate-data/test/genomes' if test else f'{SRC_DIR}/resources/genomes'
	return f'{parent_dir}/{gset}/genomes.csv'

def get_genomes_fasta_dir(wildcards_or_genomeset):
	gset = _get_genomeset(wildcards_or_genomeset)
	return f'{DL_RESOURCES}/genomes/{gset}/fasta/'

def get_genome_fasta_files(wildcards_or_genomeset, test=TEST, full_path=True):
	"""Get paths of FASTA files for the given genome set (relative to root directory)."""
	gset = _get_genomeset(wildcards_or_genomeset, require=True)

	# Ensure we've created the truncated version of this file
	if test:
		checkpoints.truncated_genome_list.get(genomeset=gset)

	parent_dir = get_genomes_fasta_dir(gset) if full_path else ''
	filenames = read_lines(get_genomes_list_file(gset, test=test))
	return [os.path.join(parent_dir, fname) for fname in filenames]


def fetch_genome_fasta_files(items, dl_dir, out_dir, nworkers, show_progress):
	"""Download FASTA files to dl_dir, then link out_dir to dl_dir."""
	from gambit_pub.download import download_parallel
	from gambit_pub.utils import symlink_to_relative

	dl_dir = Path(dl_dir)
	dl_dir.mkdir(parents=True, exist_ok=True)
	download_parallel(items, dl_dir, nworkers=nworkers, progress=show_progress)
	symlink_to_relative(dl_dir, out_dir)


# Get genomes from sets 1 and 2, from NCBI FTP server
rule fetch_genome_set_12:
	input: get_genomes_table_file
	output: directory(f'{DL_RESOURCES}/genomes/{{genomeset}}/fasta')
	params:
		dl_dir=f'{GENOMES_DL_DIR}/{{genomeset}}',
		nworkers=config['dl_nworkers'],
		show_progress=config['show_progress'],
	wildcard_constraints:
		genomeset="set[12]",
	run:
		table = pd.read_csv(input[0])
		items = [
			(NCBI_FTP_PREFIX + row.ftp_path, row.assembly_accession + '.fa.gz', row.md5)
			for _, row in table.iterrows()
		]
		fetch_genome_fasta_files(items, params['dl_dir'], output[0], params['nworkers'], params['show_progress'])


# Get genomes from sets 3 and 4, urls from config file
rule fetch_genome_set_34:
	input: get_genomes_list_file
	output: directory(f'{DL_RESOURCES}/genomes/{{genomeset}}/fasta')
	params:
		dl_dir=f'{GENOMES_DL_DIR}/{{genomeset}}',
		nworkers=config['dl_nworkers'],
		show_progress=config['show_progress'],
	wildcard_constraints:
		genomeset="set(3a|3b|4)",
	run:
		from gambit.util.io import read_lines
		files = read_lines(input[0])

		gs_dir = config['src_data']['genome_sets'][wildcards.genomeset].rstrip('/')
		prefix = GCS_PREFIX + gs_dir + '/'
		items = [(prefix + fname, fname, None) for fname in files]
		fetch_genome_fasta_files(items, params['dl_dir'], output[0], params['nworkers'], params['show_progress'])


# Download FASTQ files for genome set 3
# Unlike the FASTA download rules, this one operates on one FASTQ file at a time because they are
# large and not every one is used.
rule fetch_genome_set_3_fastq:
	output:
		f'{DL_RESOURCES}/genomes/set3/fastq/{{genome}}.fastq.gz'
	params:
		gs_dir=GCS_PREFIX + config['src_data']['genome_sets']['set3_fastq'].rstrip('/'),
	shell:
		"curl {params[gs_dir]}/{wildcards[genome]}.fastq.gz -o {output}"


def get_genome_set_3_fastq_files(wildcards=None):
	from gambit_pub.utils import stripext
	genomes = list(map(stripext, get_genome_fasta_files('set3', full_path=False)))
	return expand(rules.fetch_genome_set_3_fastq.output, genome=genomes)

# Download FASTQ files for all set 3 genomes
# Not needed by main figures and tables and not included in fetch_src_data rule
rule fetch_genome_set_3_fastq_all:
	input: get_genome_set_3_fastq_files


# Download FASTA files for genome set 5
rule fetch_genome_set_5:
	output:
		directory(f'{DL_RESOURCES}/genomes/set5/fasta'),
	params:
		url=GCS_PREFIX + config['src_data']['genome_sets']['set5'],
	shell:
		"""
		mkdir {output}
		curl {params[url]} | tar -xzf - -C {output} --strip-components=1
		"""


# Download all source data needed for main rules
rule fetch_src_data:
	input:
		*rules.fetch_gambit_db.output,
		*map(get_genomes_fasta_dir, ALL_GENOME_SETS),
		# Only the needed set 3 FASTQ files
		*expand(rules.fetch_genome_set_3_fastq.output, genome=config['figure_3']['genomes']),
