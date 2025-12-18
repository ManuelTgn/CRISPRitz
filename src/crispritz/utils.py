""" """

from pathlib import Path

import os

def validate_directory(path: str, create: bool = False) -> Path:
    dir_path = Path(path).resolve()
    if not dir_path.exists():
        if create:
            dir_path.mkdir(parents=True, exist_ok=True)
        else:
            raise FileNotFoundError(f"Directory not found: {path}")
    if not dir_path.is_dir():
        raise NotADirectoryError(f"Not a directory: {path}")
    return dir_path

def find_fasta_index(fasta_path: str) -> bool:
    return os.path.isfile(f"{fasta_path}.fai")

