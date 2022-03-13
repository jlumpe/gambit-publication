"""Utility functions to be run from the Snakemake workflow."""

import os
from pathlib import Path
from urllib.request import Request, urlopen, urlretrieve
import shutil
from tarfile import TarFile
from gzip import GzipFile
import hashlib
from concurrent.futures import ThreadPoolExecutor, as_completed

from tqdm import tqdm


def download_gcs(file: str, dest: str, untar: bool = False, gz: bool = False):
	"""Download private GCS file.

	Requires the environment variable GCS_OAUTH2_TOKEN to be set.

	Parameters
	----------
	file
		GCS file to download, {bucket}/{path}
	dest
		Destination file path.
	untar
		Use tar to extract the file with ``dest`` as the destination directory.
	gz
		If the file is gzipped an should be uncompressed.
	"""
	auth = os.environ.get('GCS_OAUTH2_TOKEN')
	if not auth:
		raise RuntimeError('Environment variable GCS_OAUTH2_TOKEN must be set to download private GCS files.')

	req = Request(f'https://storage.googleapis.com/{file}')
	req.add_header('Authorization', f'Bearer {auth}')
	fsrc = urlopen(req)
	if gz:
		fsrc = GzipFile(fileobj=fsrc)

	if untar:
		tf = TarFile(fileobj=fsrc)
		tf.extractall(dest)
	else:
		with open(dest, 'wb') as fdst:
			shutil.copyfileobj(fsrc, fdst)


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

	try:
		urlretrieve(url, file)
	except Exception as e:
		return f'Download failed: {e}'

	if checksum is not None:
		with open(file, 'rb') as f:
			md5 = get_md5(f)

		if md5 != checksum:
			file.unlink()
			return 'Download appeared to complete successfully but checksum does not match.'

	return None

def download_item(url, file, checksum=None, attempts=3):
	"""Download single file.

	Returns
	-------
	tuple
		``(success, exists, messages)`` tuple.
	"""
	file = Path(file)
	messages = []

	if file.exists():
		with open(file, 'rb') as f:
			md5 = get_md5(f)

		if checksum is not None and md5 != checksum:
			messages += ['File exists but checksum is incorrect (partial download?). Deleting.']
			file.unlink()
		else:
			return True, True, messages

	for i in range(attempts):
		last_err = attempt_download(url, file, checksum)
		if not last_err:
			return True, False, messages

	messages.append(f'Download failed after {attempts} attempts. Last error: {last_err}')
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

	with ThreadPoolExecutor(max_workers=nworkers) as executor:
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
