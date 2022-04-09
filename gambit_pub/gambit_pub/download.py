"""Functions to download and validate source files."""

import os
from pathlib import Path
from urllib.request import urlretrieve
from urllib.error import HTTPError
import hashlib
from concurrent.futures import ProcessPoolExecutor, as_completed
from time import sleep

from tqdm import tqdm


def get_md5(file):
	"""Calculate MD5 hash of a readable file object."""
	md5 = hashlib.md5()

	while True:
		chunk = file.read(8192)
		if not chunk:
			break
		md5.update(chunk)

	return md5.hexdigest()


def attempt_download(url, file, checksum=None):
	"""Attempt to download URL to file.

	Returns
	-------
	(err_msg, should_retry) tuple.
	"""

	try:
		urlretrieve(url, file)
	except Exception as e:
		if file.is_file():
			file.unlink()

		should_retry = not (isinstance(e, HTTPError) and 400 <= e.code < 500)
		return f'Download failed: {e}', should_retry

	if checksum is not None:
		with open(file, 'rb') as f:
			md5 = get_md5(f)

		if md5 != checksum:
			file.unlink()
			return 'Download appeared to complete successfully but checksum does not match.', True

	return None, False

def download_item(url, file, checksum=None, attempts=3, retry_delay=.25):
	"""Download single file.

	Returns
	-------
	tuple
		``(success, exists, messages)`` tuple.
	"""
	file = Path(file)
	messages = []

	if file.exists():

		if checksum is not None:
			with open(file, 'rb') as f:
				md5 = get_md5(f)

			if md5 != checksum:
				messages += ['File exists but checksum is incorrect (partial download?). Deleting.']
				file.unlink()
			else:
				return True, True, messages
		else:
			return True, True, messages

	for i in range(attempts):
		if i > 0:
			sleep(retry_delay)

		last_err, should_retry = attempt_download(url, file, checksum)

		# Success
		if not last_err:
			return True, False, messages

		if not should_retry:
			messages.append(last_err)
			break

	else:
		# Went through all attempts, no success
		messages.append(f'Failed after {attempts} attempts. Last error: {last_err}')

	return False, False, messages


def download_parallel(items, out_dir, nworkers=None):
	"""Download multiple files in parallel.

	Parameters
	----------
	items:
		``(url, file, checksum)`` tuples.
	out_dir
	nworkers: int
	"""
	out_dir = Path(out_dir)
	if nworkers is None:
		nworkers = 2 * os.cpu_count()

	nsuccess = 0
	nexists = 0
	nfailed = 0

	with ProcessPoolExecutor(max_workers=nworkers) as executor:
		future_to_file = {
			executor.submit(download_item, url, out_dir / file, checksum): file
			for url, file, checksum in items
		}

		with tqdm(as_completed(future_to_file), total=len(future_to_file)) as pbar:
			for future in pbar:
				file = future_to_file[future]
				success, exists, messages = future.result()

				for msg in messages:
					pbar.write(f'{file}: {msg}')

				if exists:
					nexists += 1
				elif success:
					nsuccess += 1
				else:
					nfailed += 1

	# return nsuccess, nexists, nfailed
	print(f'{nsuccess} files successfully downloaded, {nexists} already exist, {nfailed} failed.')
	if nfailed > 0:
		raise RuntimeError(f'{nfailed} downloads failed.')
