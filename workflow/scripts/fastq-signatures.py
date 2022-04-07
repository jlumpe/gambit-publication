"""
Generate k-mer signatures from a FASTQ file with a range of parameter values, and calculate Jaccard
distances to the fasta_sig derived from the assembled FASTA file.

This uses custom code to do all the parameters at once, so we're not using the GAMBIT CLI or
standard API functions for fasta_sig and distance calculation.

Expected Snakemake variables:
* input
  * signatures: Signature file derived from FASTA files.
  * fastq: FASTQ file.
* params
  * fasta_name: Name of FASTA file to look up the correct signature.
* output: CSV file to write results to.
"""

import numpy as np
import pandas as pd

from gambit.sigs import load_signatures
from gambit.seq import SequenceFile
from gambit_pub.fastq import PhredAccumulator, accumulate_kmers_fastq, phredsum


### Parameters ###

# Way too many bins for plot, but doesn't cost much and we can re-bin later.
MIN_PHRED = np.arange(0, 41)
MIN_COUNT = np.arange(100) + 1

AGG_FUNCS = [np.min, phredsum]
AGG_NAMES = ['min', 'phredsum']


### Get signatures of FASTA file ###

fasta_sigs = load_signatures(snakemake.input['signatures'])
kmerspec = fasta_sigs.kmerspec

fasta_filename = snakemake.params['fasta_name']
(fasta_sig_index,) = np.flatnonzero(fasta_sigs.ids == fasta_filename)
fasta_sig = fasta_sigs[fasta_sig_index]
fasta_sig_len = len(fasta_sig)


### Count kmers in fastq file ###

seqfile = SequenceFile(snakemake.input['fastq'], 'fastq', 'gzip')

phred_bins = MIN_PHRED[1:]  # Left bin edges, count anything <1 as 0
accums = [PhredAccumulator(phred_bins) for _ in AGG_FUNCS]

for record in seqfile.parse():
    accumulate_kmers_fastq(kmerspec, record, zip(AGG_FUNCS, accums))


### Statistics ###

phred_index = pd.Series(MIN_PHRED, name='min_phred')
count_index = pd.Series(MIN_COUNT, name='min_count')

out_chunks = []

for agg_name, accum in zip(AGG_NAMES, accums):

    # Indices of found k-mers, and the number found in each PHRED bin
    kmer_inds, kmer_counts = accum.to_arrays()

    # Number of k-mers common to fasta and fastq signatures, or just in fastq
    # by minimum count and minimum aggregated PHRED score
    out_shape = (len(MIN_COUNT), len(MIN_PHRED))
    n_intersect = np.zeros(out_shape, dtype=np.int32)
    n_fastq = np.zeros(out_shape, dtype=np.int32)

    # For each k-mer found
    for i, idx in enumerate(kmer_inds):
        # Which array to increment
        a = n_intersect if idx in fasta_sig else n_fastq

        counts_cum = np.cumsum(kmer_counts[i, ::-1])[::-1]
        for j, min_count in enumerate(MIN_COUNT):
            a[j, :] += counts_cum >= min_count

    # To tall-skinny table format
    chunk = pd.concat(
        [
            pd.DataFrame(arr, index=count_index, columns=phred_index).unstack()
            for arr in [n_intersect, n_fastq]
        ],
        axis=1,
        keys=['intersection', 'fastq_only'],
    )

    chunk['fasta_only'] = fasta_sig_len - chunk['intersection']
    n_union = chunk['fastq_only'] + fasta_sig_len
    chunk['jaccard'] = (chunk['fasta_only'] + chunk['fastq_only']) / n_union

    out_chunks.append(chunk)


### Write output ###

out_df = pd.concat(out_chunks, keys=AGG_NAMES, names=['aggregation'])
out_df.to_csv(snakemake.output[0])
