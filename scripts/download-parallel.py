#!/usr/bin/env python3

"""
Download multiple files in parallel, with the list of URLs taken from a .csv file, and optionally
verify MD5 checksums.

Mostly just used to download genomes from the NCBI FTP server.
"""


import os
import sys
import argparse
from urllib.request import urlretrieve
from pathlib import Path
import hashlib
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
from tqdm import tqdm


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--url', default='url', help='Column of CSV file containing the URL to download')
parser.add_argument('--file', help='Column of CSV file containing the destination file name')
parser.add_argument('--md5', help='Column of CSV file containing expected MD5 checksum')
parser.add_argument('--id', help='Column of CSV file containing unique ID for error reporting')
parser.add_argument('-o', '--out-dir', default='.', help='Directory to write files to')
parser.add_argument('-t', '--threads', help='Number of threads to use')
parser.add_argument('csvfile', metavar='CSVFILE', help='CSV file containing URLs')


def main():
	args = parser.parse_args()

	df = pd.read_csv(args.csvfile)
	items = make_items(df, args.url, args.file, args.md5, args.id, args.out_dir)

	nsuccess, nexists, nfailed = download_items(items, args.threads)

	print(f'{nsuccess} files successfully downloaded, {nexists} already exist, {nfailed} failed.')
	sys.exit(1 if nfailed > 0 else 0)


def make_items(df, url_col, file_col=None, md5_col=None, id_col=None, out_dir='.'):
	items = []

	for i, row in df.iterrows():
		url = row[url_col]
		id_ = url if id_col is None else row[id_col]
		file = Path(out_dir) / (url.rsplit('/', 1)[1] if file_col is None else row[file_col])
		md5 = None if md5_col is None else row[md5_col]
		items.append((id_, url, file, md5))

	return items


def get_md5(file):
	"""Calculate MD5 hash of an open file object."""
	md5 = hashlib.md5()

	while True:
		chunk = file.read(8192)
		if not chunk:
			break
		md5.update(chunk)

	return md5.hexdigest()


def download_item(url, file, checksum=None):
	"""Download single file.

	Returns
	-------
	tuple
		``(success, exists, messages)`` tuple.
	"""

	messages = []

	if file.exists():
		with open(file, 'rb') as f:
			md5 = get_md5(f)

		if checksum is not None and md5 != checksum:
			messages += ['File exists but checksum is incorrect (partial download?). Deleting.']
			file.unlink()
		else:
			return True, True, messages

	try:
		urlretrieve(url, file)
	except Exception as e:
		messages.append(f'Download failed: {e}')
		return False, False, messages

	if checksum is not None:
		with open(file, 'rb') as f:
			md5 = get_md5(f)

		if md5 != checksum:
			messages.append('Download appeared to complete successfully but checksum does not match.')
			file.unlink()
			return False, False, messages

	return True, False, messages


def download_items(items, nworkers=None):
	"""Download multiple files.

	Parameters
	----------
	items:
		``(id, url, file, checksum)`` tuples.
	nworkers: int
	"""

	if nworkers is None:
		nworkers = 2 * os.cpu_count()

	nsuccess = 0
	nexists = 0
	nfailed = 0

	with ThreadPoolExecutor(max_workers=nworkers) as executor:
		future_to_id = {
			executor.submit(download_item, url, file, checksum): id_
			for id_, url, file, checksum in items
		}

		with tqdm(as_completed(future_to_id), total=len(future_to_id)) as pbar:
			for future in pbar:
				id_ = future_to_id[future]
				success, exists, messages = future.result()

				for msg in messages:
					pbar.write(f'{id_} {msg}')

				if success:
					if exists:
						nexists += 1
					else:
						nsuccess += 1
				else:
					nfailed += 1

	return nsuccess, nexists, nfailed


if __name__ == '__main__':
	main()
