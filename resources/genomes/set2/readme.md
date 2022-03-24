# Genome set 2

Set of 70 genomes derived from Supporting Table 1 in (Konstantinidis 2016). This is used as a
high-diversity data set for evaluating the GAMBIT distance metric as it spans multiple bacterial
phyla. We manually selected the best RefSeq assembly to represent each strain.


## Files

* `genomes.csv` - Table of metadata for all genomes.
  * `group` - Value of "group" column in original table.
  * `strain` - Value of "strain" column in original table.
  * `assembly_uid` - UID of selected NCBI assembly database entry.
  * `assembly_accession` - Accession number of selected NCBI assembly database entry.
  * `organism` - From NCBI database entry metadata.
  * `taxid` - From NCBI database entry metadata.
  * `filename` - Name of sequence file in `fasta/` directory.
  * `url` - URL of sequence file on NCBI FTP server.
  * `md5` - MD5 checksum of sequence file.
* `fasta/` - Genome FASTA files are downloaded to this directory.


## References

Konstantinidis, K. T., & Tiedje, J. M. (2005). Genomic insights that advance the species definition for prokaryotes. Proceedings of the National Academy of Sciences of the United States of America, 102(7), 2567–2572.
