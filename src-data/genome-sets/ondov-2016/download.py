# Downloads sequences from NCBI FTP server to fata/ directory.
# Run from this directory, not the repo's root.

import sys
import os
from urllib.request import urlretrieve
from pathlib import Path
import hashlib
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
from tqdm import tqdm


def get_md5(file):
	"""Calculate MD5 hash of an open file object."""
	md5 = hashlib.md5()

	while True:
		chunk = file.read(8192)
		if not chunk:
			break
		md5.update(chunk)

	return md5.hexdigest()


def download_seq(url, file, md5hash):
	messages = []

	if file.exists():
		with open(file, 'rb') as f:
			file_hash = get_md5(f)

		if file_hash != md5hash:
			messages += ['File exists but hash is incorrect (partial download?). Deleting.']
			file.unlink()
		else:
			return True, True, messages

	try:
		urlretrieve(url, file)
	except Exception as e:
		messages.append(f'Download failed: {e}')
		return False, False, messages

	with open(file, 'rb') as f:
		file_hash = get_md5(f)

	if file_hash != md5hash:
		messages.append('Download appeared to complete successfully but file hash does not match.')
		# file.unlink()
		return False, False, messages

	return True, False, messages


if __name__ == '__main__':
	genomes = pd.read_csv('genomes.csv')
	genomes['file'] = list(map(Path, genomes['file']))

	nsuccess = 0
	nexists = 0
	nfailed = 0

	with ThreadPoolExecutor(max_workers=8) as executor:
		future_to_acc = {
			executor.submit(download_seq, row.url, row.file, row.md5): row.assembly_accession
			for idx, row in genomes.iterrows()
		}

		with tqdm(as_completed(future_to_acc), total=len(future_to_acc)) as pbar:
			for future in pbar:
				acc = future_to_acc[future]
				success, exists, messages = future.result()

				for msg in messages:
					pbar.write(f'{acc} {msg}')

				if success:
					if exists:
						nexists += 1
					else:
						nsuccess += 1
				else:
					nfailed += 1

	print(f'{nsuccess} genomes successfully downloaded, {nexists} already exist, {nfailed} failed.')
	sys.exit(1 if nfailed > 0 else 0)
