# Download source data


# Download GAMBIT database files
rule get_gambit_db:
	output:
		genomes='resources/gambit-db/db-genomes.db',
		signatures='resources/gambit-db/db-signatures.h5',
	shell:
	    '''
		wget -O {output.genomes} {config[gambit_db][base_url]}/{config[gambit_db][genomes_file]}
		wget -O {output.signatures} {config[gambit_db][base_url]}/{config[gambit_db][signatures_file]}
		'''


# Download ondov-2016 genomes
rule get_genomes_ondov2016:
	output:
		directory('resources/genomes/ondov_2016/fasta'),
	run:
		os.mkdir(output[0])
		table = pd.read_csv('resources/genomes/ondov_2016/genomes.csv')
		items = [(row.url, row.assembly_accession + '.fa.gz', row.md5) for _, row in table.iterrows()]
		wf_utils.download_parallel(items, output[0])


# Download konstantinidis-2005 genomes
rule get_genomes_konstantinidis2005:
	output:
		directory('resources/genomes/konstantinidis_2005/fasta'),
	run:
		os.mkdir(output[0])
		table = pd.read_csv('resources/genomes/konstantinidis_2005/genomes.csv')
		items = [(row.url, row.assembly_accession + '.fa.gz', row.md5) for _, row in table.iterrows()]
		wf_utils.download_parallel(items, output[0])


# Download genomes for figure 6
rule get_genomes_fig6:
	output:
	    directory('resources/genomes/figure_6/fasta'),
	run:
		wf_utils.download_gcs(
			'hesslab-gambit/genomes/210910-ecoli-genomes-for-tree/fasta.tar.gz',
			Path(output[0]).parent,
			untar=True,
			gz=True,
		)


# Download all source data
rule get_src_data:
	input:
	     rules.get_gambit_db.output,
	     rules.get_genomes_konstantinidis2005.output,
	     rules.get_genomes_ondov2016.output,
	     rules.get_genomes_fig6.output,
