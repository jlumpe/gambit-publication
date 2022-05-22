"""Get k-mer signatures from FASTQ files using a a range of parameters."""

import numpy as np

from gambit.kmers import find_kmers


def get_phred(record):
	return np.asarray(record.letter_annotations['phred_quality'])


class PhredAccumulator:
	"""Accumulate k-mer counts binned by aggregated PHRED score."""

	def __init__(self, bins_left, dtype=np.dtype('u2')):
		self.bins_left = np.asarray(bins_left)
		assert np.all(np.diff(self.bins_left) > 0), 'Bins must be increasing'
		self.nbins = len(self.bins_left)
		self.dtype = dtype
		self.dict = dict()

	def get_bin(self, score):
		return np.searchsorted(self.bins_left, score, side='right') - 1

	def add(self, index, score):
		"""Add k-mer by index and minimum PHRED score."""
		b = self.get_bin(score)
		if b < 0:  # Less than smallest bin
			return

		try:
			arr = self.dict[index]
		except KeyError:
			arr = self.dict[index] = np.zeros(self.nbins, dtype=self.dtype)

		arr[b] += 1
		assert arr[b] != 0  # catch overflow

	def to_arrays(self):
		"""Convert to array of found k-mer indices and array of counts."""
		indices = np.fromiter(self.dict, dtype=int)
		indices.sort()

		counts = np.empty((len(indices), self.nbins), dtype=self.dtype)
		for row, index in enumerate(indices):
			counts[row, :] = self.dict[index]

		return indices, counts


def accumulate_kmers_fastq(kspec, record, accumulator):
	"""Add k-mers in single FASTQ record to accumulator."""
	phred = get_phred(record)

	for match in find_kmers(kspec, record.seq):
		try:
			index = match.kmer_index()
		except ValueError:
			continue  # Invalid nucleotide code, ignore

		min_phred = np.min(phred[match.full_indices()])
		accumulator.add(index, min_phred)
