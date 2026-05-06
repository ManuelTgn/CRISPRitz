""" """

from pathlib import Path

import sys
import os

# ------------------------------------------------------------------------------
# define static variables
# ------------------------------------------------------------------------------

# toolname and command
TOOLNAME = "CRISPRitz"
COMMAND = "crispritz"

# sub commands
SUBCOMMANDS = ["add-variants", "index-genome", "search"]


# ------------------------------------------------------------------------------
# define utility functions
# ------------------------------------------------------------------------------


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


def create_folder(folder: str) -> str:
    os.makedirs(folder, exist_ok=True)
    assert os.path.isdir(folder)
    return os.path.abspath(folder)


def find_fasta_index(fasta_path: str) -> bool:
    return os.path.isfile(f"{fasta_path}.fai")


def find_tabix_index(fname_path: str) -> bool:
    return os.path.isfile(f"{fname_path}.tbi")


def set_processes(tasks_n: int, threads: int) -> int:
    return min(tasks_n, threads)
