"""
Generate supplemental figures and data.
"""


rule sfig1:
	input: expand(rules.gambit_ani_correlation.output, paramspace='full')
	output: 'results/figures/supplemental-figure-1.png'
	params:
		paramspace='full',
	script: '../scripts/figure-2.py'


rule sfig2:
	input: expand(rules.gambit_vs_ani.output['pairs'], k=K, prefix=PREFIX, genomeset=COMPARISON_GENOME_SETS)
	output:
		figure='results/figures/supplemental-figure-2.png',
	params:
		genome_sets=COMPARISON_GENOME_SETS,
	script: '../scripts/supplemental-figure-2.py'


# Supplemental table 3: reference database genomes
rule stable3:
	input: rules.fetch_gambit_db.output['genomes']
	output: 'results/tables/supplemental-table-3.csv'
	run:
		import pandas as pd
		from gambit.db import load_genomeset
		from gambit_pub.utils import getattr_coalesce, fix_int_cols

		session, gset = load_genomeset(input[0])
		rows = []

		with session:
			for genome in gset.genomes:
				species = genome.taxon.ancestor_of_rank('species')
				genus = genome.taxon.ancestor_of_rank('genus')
				rows.append(dict(
					db_id=genome.genome_id,
					assembly_accession=genome.refseq_acc,
					assembly_uid=genome.ncbi_id,
					genus_name=getattr_coalesce(genus, 'name'),
					genus_db_id=getattr_coalesce(genus, 'id'),
					genus_ncbi_id=getattr_coalesce(genus, 'ncbi_id'),
					species_name=species.name,
					species_db_id=species.id,
					species_ncbi_id=species.ncbi_id,
				))

		df = pd.DataFrame(rows)
		fix_int_cols(df, ['genus_db_id', 'genus_ncbi_id', 'species_db_id', 'species_ncbi_id'])
		df.sort_values(['genus_ncbi_id', 'species_ncbi_id', 'db_id'], inplace=True)
		df.to_csv(output[0], index=False)


# Supplemental table 4: reference database taxa
rule stable4:
	input: rules.fetch_gambit_db.output['genomes']
	output: 'results/tables/supplemental-table-4.csv'
	run:
		import pandas as pd
		from gambit.db import load_genomeset
		from gambit_pub.utils import fix_int_cols

		session, gset = load_genomeset(input[0])
		rows = []

		with session:
			for taxon in gset.taxa:
				if taxon.rank is not None:
					genus = taxon.ancestor_of_rank('genus')
					rows.append(dict(
						db_id=taxon.id,
						name=taxon.name,
						rank=taxon.rank,
						parent_id=genus.id if genus not in (taxon, None) else None,
						ncbi_id=taxon.ncbi_id,
						threshold=taxon.distance_threshold,
					))

		df = pd.DataFrame(rows)
		fix_int_cols(df, ['parent_id', 'ncbi_id'])
		df.sort_values(['rank', 'name'], inplace=True)
		df.to_csv(output[0], index=False)
