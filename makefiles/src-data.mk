# Fetch source data


.PHONY: \
  src-data \
  clean-src-data \
  src-data/gambit/database \
  src-data/genomes/ondov-2016/fasta \
  clean--src-data/gambit/database \
  clean--src-data/genomes/ondov-2016/fasta
.ONESHELL: \
  src-data/genomes/figure-6/fasta


# Get all source data
src-data: \
  src-data/gambit/database \
  src-data/genomes/ondov-2016/fasta

# Remove all source data
clean-src-data: \
  clean--src-data/gambit/database \
  clean--src-data/genomes/ondov-2016/fasta


# GAMBIT database files
src-data/gambit/database: src-data/gambit/database/db-genomes.db src-data/gambit/database/db-signatures.h5

_gambit_db_base_url := https://storage.googleapis.com/hesslab-gambit-public/databases/refseq-curated/1.0-beta

src-data/gambit/database/db-genomes.db:
	# Fetch GAMBIT database genomes
	wget -O $@ $(_gambit_db_base_url)/gambit-genomes-1.0b1-210719.db

src-data/gambit/database/db-signatures.h5:
	# Fetch GAMBIT database signatures
	wget -O $@ $(_gambit_db_base_url)/gambit-signatures-1.0b1-210719.h5

clean--src-data/gambit/database:
	rm src-data/gambit/database/db-genomes.db
	rm src-data/gambit/database/db-signatures.h5


# Ondov 2016 genomes
src-data/genomes/ondov-2016/fasta: src-data/genomes/ondov-2016/fasta/.completed

src-data/genomes/ondov-2016/fasta/.completed:
	# Download Ondov 2016 genomes
	mkdir -p $(dir $@)
	$(conda run) python scripts/download-parallel.py \
		--id=assembly_accession --file=file --md5=md5 \
		-o src-data/genomes/ondov-2016/fasta \
		src-data/genomes/ondov-2016/genomes.csv
	touch $@

clean--src-data/genomes/ondov-2016/fasta:
	rm -r src-data/genomes/ondov-2016/fasta


# Figure 6 genomes
src-data/genomes/figure-6/fasta:
	# Download figure 6 genomes
	cd $(dir $@)
	$(call DOWNLOAD_GCS,hesslab-gambit/genomes/210910-ecoli-genomes-for-tree/fasta.tar.gz,kjdkfj)
	tar -xzf fasta.tar.gz
	rm fasta.tar.gz
