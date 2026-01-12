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

# verbosity levels
VERBOSITYLVL = [0, 1, 2, 3]

# sub commands
SUBCOMMANDS = ["add-variants"]


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


def print_verbosity(message: str, verbosity: int, verbosity_threshold: int) -> None:
    """Print a message if the verbosity level meets the threshold.

    Writes the message to standard output if the current verbosity is greater
    than or equal to the specified threshold.

    Args:
        message (str): The message to print.
        verbosity (int): The current verbosity level.
        verbosity_threshold (int): The minimum verbosity level required to print
            the message.

    Returns:
        None
    """
    if verbosity >= verbosity_threshold:
        sys.stdout.write(f"{message}\n")
    return


def create_folder(folder: str) -> str:
    os.makedirs(folder, exist_ok=True)
    assert os.path.isdir(folder)
    return os.path.abspath(folder)


def find_fasta_index(fasta_path: str) -> bool:
    return os.path.isfile(f"{fasta_path}.fai")
