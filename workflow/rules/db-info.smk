"""
Rules which derive some simple information and statistics from the GAMBIT reference database.
"""

rule db_basic_stats:
	input:
		genomes=rules.fetch_gambit_db.output['genomes'],
		signatures=rules.fetch_gambit_db.output['signatures'],
	output: 'results/db-info/basic-stats.txt'
	run:
		from gambit.db import ReferenceDatabase
		db = ReferenceDatabase.load(input['genomes'], input['signatures'])
		stats = dict()

		stats['Genome count'] = len(db.genomes)
		stats['Species count'] = db.genomeset.taxa.filter_by(rank='species').count()
		stats['Genus count'] = db.genomeset.taxa.filter_by(rank='genus').count()

		# Prevents SQLA printing errors due to whatever the hell Snakemake is doing with threads.
		db.session.close()

		with open(output[0], 'w') as f:
			for label, value in stats.items():
				print(label, value, sep='\t', file=f)


rule db_kmer_counts:
	input:
		signatures=rules.fetch_gambit_db.output['signatures'],
	output: 'results/db-info/kmer-counts-summary.tsv'
	run:
		import pandas as pd
		from gambit.sigs import load_signatures
		sigs = load_signatures(input['signatures'])
		summary = pd.Series(sigs.sizes()).describe()
		summary.to_csv(output[0], header=False)


rule db_info:
	input:
		*rules.db_basic_stats.output,
		*rules.db_kmer_counts.output,
