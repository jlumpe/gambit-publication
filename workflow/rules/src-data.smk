"""Download source data."""


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


# Download genome set 1 or 2
rule get_genome_set_12:
	output:
		directory("resources/genomes/{genomeset}/fasta/"),
	wildcard_constraints:
		genomeset="set[12]",
	run:
		from gambit_pub.download import download_parallel

		outdir = Path(output[0])
		dl_dir = outdir.parent / '.fasta-download'
		dl_dir.mkdir(exist_ok=True)

		table = pd.read_csv(outdir.parent / 'genomes.csv')
		items = [(row.url, row.assembly_accession + '.fa.gz', row.md5) for _, row in table.iterrows()]
		download_parallel(items, dl_dir)

		outdir.symlink_to(dl_dir.name, True)


# Download genome set 3
rule get_genome_set_3:
	output:
		directory("resources/genomes/set3/fasta/")
	run:
		from gambit_pub.download import download_parallel

		gs_dir = config['src_data']['genome_sets']['set3']['fasta'].rstrip('/')
		prefix = GCS_PREFIX + gs_dir + '/'

		outdir = Path(output[0])
		dl_dir = outdir.parent / '.fasta-download'
		dl_dir.mkdir(exist_ok=True)

		items = []
		with open('resources/genomes/set3/genomes.txt') as f:
			for line in f:
				fname = line.strip()
				items.append((prefix + fname, fname, None))

		download_parallel(items, dl_dir)
		outdir.symlink_to(dl_dir.name, True)


# Download fastq files for genome set 3
rule get_genome_set_3_fastq:
	output:
	      directory("resources/genomes/set3/fastq/")
	run:
		from gambit_pub.download import download_parallel

		gs_dir = config['src_data']['genome_sets']['set3']['fastq'].rstrip('/')
		prefix = GCS_PREFIX + gs_dir + '/'

		outdir = Path(output[0])
		dl_dir = outdir.parent / '.fasta-download'
		dl_dir.mkdir(exist_ok=True)

		items = []
		with open('resources/genomes/set3/genomes.txt') as f:
			for line in f:
				fname = line.strip().rsplit('.', 1)[0] + '.fasta.gz'
				items.append((prefix + fname, fname, None))

		download_parallel(items, dl_dir)
		outdir.symlink_to(dl_dir.name, True)


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
	     *expand('resources/genomes/{gset}/fasta', gset=['set1', 'set2', 'set3', 'figure_6']),
	     rules.get_genome_set_3_fastq.output,
