# Ondov 2016 Escherichia Genomes

Set of 492 Escherichia genomes used in the Mash paper (Ondov 2016). This is used as a
medium-diversity data set for evaluating the GAMBIT distance metric.

Genomes are RefSeq assemblies from NCBI. List is available from the [documentation
page](https://mash.readthedocs.io/en/latest/data.html) in the `Escherichia.tar.gz` archive. We were
unable to determine the correct accession number for 8 of the genomes in this list, which is why
there are 492 genomes in this set rather than 500.


## Files

* `genomes.csv` - Table of metadata for all genomes.
* `fasta/` - Genome FASTA files are downloaded to this directory.


## References

Ondov, B. D., Treangen, T. J., Melsted, P., Mallonee, A. B., Bergman, N. H., Koren, S., & Phillippy, A. M. (2016). Mash: fast genome and metagenome distance estimation using MinHash. Genome Biology, 17(1), 132.
