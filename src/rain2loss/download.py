"""Download helper: skips files that already exist and shows the progress in MB."""
import urllib.request
from pathlib import Path


def download(url, path):
    """Download url to path unless the file already exists."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        print("already there:", path.name)
        return

    def progress(blocks, block_size, total):
        print(f"\r{blocks * block_size / 1e6:8.1f} / {total / 1e6:.1f} MB", end="")

    partial = path.with_name(path.name + ".part")
    urllib.request.urlretrieve(url, partial, reporthook=progress)
    partial.replace(path)                       # rename only when the download is complete
    print("\ndone:", path.name)