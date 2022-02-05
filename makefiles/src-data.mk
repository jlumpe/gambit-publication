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
#  src-data/genomes/ondov-2016/fasta \
# Don't delete data if the rule fails, probably just want to restart the download
.PRECIOUS: src-data/genomes/ondov-2016/fasta


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

src-data/gambit/db-genomes.db:
	# Fetch GAMBIT database genomes
	wget -O $@ $(GAMBIT_DB_GENOMES_URL)

src-data/gambit/db-signatures.h5:
	# Fetch GAMBIT database signatures
	wget -O $@ $(GAMBIT_DB_SIGNATURES_URL)

clean--src-data/gambit/database:
	rm src-data/gambit/database/db-genomes.db
	rm src-data/gambit/database/db-signatures.h5


# Ondov 2016 genomes
src-data/genomes/ondov-2016/fasta: src-data/genomes/ondov-2016/fasta/.completed

src-data/genomes/ondov-2016/fasta/.completed:
	# Download Ondov 2016 genomes
	mkdir -p $(dir $@)
	(cd src-data/genomes/ondov-2016; $(CONDA_RUN) python download.py)
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
