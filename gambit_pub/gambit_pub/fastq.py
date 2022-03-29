"""Get k-mer signatures from FASTQ files using a a range of parameters."""

import numpy as np

from gambit.kmers import find_kmers


# Factor to convert exponent from PHRED (base 10 ** (1/10)) to natural (base e)
PHRED_TO_NAT = -np.log(10) / 10

def phredsum(q):
	return np.logaddexp.reduce(np.asarray(q) * PHRED_TO_NAT) / PHRED_TO_NAT

def get_phred(record):
	return np.asarray(record.letter_annotations['phred_quality'])

class PhredAccumulator:
	"""Accumulate k-mer counts binned by aggregated PHRED score."""

	def __init__(self, bin_edges, dtype=np.dtype('u2')):
		self.bin_edges = np.asarray(bin_edges)
		self.nbins = len(self.bin_edges) + 1
		self.dtype = dtype
		self.dict = dict()

	def get_bin(self, score):
		return np.searchsorted(self.bin_edges, score, side='right')

	def add(self, index, score):
		b = self.get_bin(score)

		try:
			arr = self.dict[index]
		except KeyError:
			arr = self.dict[index] = np.zeros(self.nbins, dtype=self.dtype)

		arr[b] += 1
		assert arr[b] != 0  # catch overflow

	def to_arrays(self):
		indices = np.fromiter(self.dict, dtype=int)
		indices.sort()

		counts = np.empty((len(indices), self.nbins), dtype=self.dtype)
		for row, index in enumerate(indices):
			counts[row, :] = self.dict[index]

		return indices, counts


def accumulate_kmers_fastq(kspec, record, accumulators):
	phred = get_phred(record)

	for match in find_kmers(kspec, record.seq):
		try:
			index = match.kmer_index()
		except ValueError:
			continue

		p = phred[match.full_indices()]

		for agg_func, accum in accumulators:
			accum.add(index, agg_func(p))
