# Fetch source data

#include src-data/gambit/targets.mk
#include src-data/genome-sets/figure-6/targets.mk


.PHONY: src-data clean-src-data src-data/gambit/database
.ONESHELL: src-data/genome-sets/figure-6/fasta

# Get all source data
src-data: src-gambit-db

# Remove all source data
clean-src-data: clean-src-gambit-db


# GAMBIT database files
src-data/gambit/database: src-data/gambit/database/db-genomes.db src-data/gambit/database/db-signatures.h5

src-data/gambit/db-genomes.db:
	# Fetch GAMBIT database genomes
	wget -O $@ $(GAMBIT_DB_GENOMES_URL)

src-data/gambit/db-signatures.h5:
	# Fetch GAMBIT database signatures
	wget -O $@ $(GAMBIT_DB_SIGNATURES_URL)


# Figure 6 data
src-data/genome-sets/figure-6/fasta:
	# Download figure 6 genomes
	cd $(dir $@)
	$(call DOWNLOAD_GCS,hesslab-gambit/genomes/210910-ecoli-genomes-for-tree/fasta.tar.gz,kjdkfj)
	tar -xzf fasta.tar.gz
	rm fasta.tar.gz
