"""
Count occurrences of k-mers in a FASTQ filtered by PHRED quality scores.

Expected Snakemake variables:
* input
  * signatures: Signature file derived from FASTA files.
  * fastq: FASTQ file.
  * genomes_table: CSV file of genome/assembly data.
* params
  * fasta_id: ID used to look up the correct FASTA signature.
  * min_phred: Minimum PHRED scores to filter reads by.
  * truncate_reads: Only process this many reads (speed up this rule for testing).
* output: Directory with the following files:
  * kmer-counts.csv: CSV file of counts for each k-mer found in FASTQ file.
  * stats.json: Additional statistics.
"""

from pathlib import Path
import json

import numpy as np
import pandas as pd

from gambit.sigs import load_signatures
from gambit.seq import SequenceFile
from gambit_pub.fastq import PhredAccumulator, accumulate_kmers_fastq


TRUNCATE_READS = snakemake.params['truncate_reads']

### Get FASTA file stats ###

fasta_df = pd.read_csv(snakemake.input['genomes_table'], index_col=0)
fasta_stats = fasta_df.loc[snakemake.wildcards['genome']]


### Get signature of FASTA file ###

fasta_sigs = load_signatures(snakemake.input['signatures'])
kmerspec = fasta_sigs.kmerspec

fasta_id = snakemake.params['fasta_id']
(fasta_sig_index,) = np.flatnonzero(fasta_sigs.ids == fasta_id)
fasta_sig = fasta_sigs[fasta_sig_index]
fasta_sig_len = len(fasta_sig)


### Count kmers in fastq file ###

seqfile = SequenceFile(snakemake.input['fastq'], 'fastq', 'gzip')

phred_bins = snakemake.params['min_phred']
accumulator = PhredAccumulator(phred_bins)

nreads = 0
total_len = 0

for record in seqfile.parse():
	nreads += 1
	total_len += len(record.seq)
	accumulate_kmers_fastq(kmerspec, record, accumulator)

	if TRUNCATE_READS is not None and nreads >= TRUNCATE_READS:
		break


### Counts table ###

indices, counts_binned = accumulator.to_arrays()
# Make counts cumulative instead of binned
counts_cum = np.cumsum(counts_binned[:, ::-1], axis=1)[:, ::-1]

counts_df = pd.DataFrame(
	counts_cum,
	index=pd.Series(indices, name='index'),
	columns=phred_bins,
)

counts_df.insert(0, 'in_fasta', np.in1d(indices, fasta_sig))


### Additional statistics ###

assembly_len = int(fasta_stats.loc['total_length'])

stats = dict(
	nreads=nreads,
	total_length=total_len,
	estimated_coverage=total_len / assembly_len,
	assembly_length=assembly_len,
	assembly_nkmers=len(fasta_sig),
)


### Write output ###

outdir = Path(snakemake.output[0])
outdir.mkdir(parents=True, exist_ok=True)

counts_df.to_csv(outdir / 'kmer-counts.csv')

with open(outdir / 'stats.json', 'w') as f:
	json.dump(stats, f)
