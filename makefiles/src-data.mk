# Fetch source data

.PHONY: src-data src-gambit-db


# Get all source data
src-data: src-gambit-db


# GAMBIT database files
src-gambit-db: src-data/gambit/db-genomes.db src-data/gambit/db-signatures.h5

src-data/gambit/db-genomes.db:
	wget -O $@ $(GAMBIT_DB_GENOMES_URL)

src-data/gambit/db-signatures.h5:
	wget -O $@ $(GAMBIT_DB_SIGNATURES_URL)
