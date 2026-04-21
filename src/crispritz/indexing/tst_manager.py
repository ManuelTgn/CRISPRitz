""" """

from ..crispritz_error import CrispritzTstError
from ..exception_handlers import exception_handler
from ..crispritz_cpp import build_tree_cpp
from ..utils import print_verbosity, VERBOSITYLVL
from ..genome_io import GenomeReader
from ..pam import PAM

from typing import List

import os


def build_ternary_search_tree(
    fastas: List[str],
    pam_file: str,
    bmax: int,
    threads: int,
    verbosity: int,
    debug: bool,
) -> None:
    # read pam file for constructing ternary search tree for current pam
    pam = PAM(pam_file, debug)
    # construct ternary search tree for each input fasta file
    for fasta in fastas:
        reader = GenomeReader(fasta, debug)  # read input enriched fasta
        reader.read()
        try:  # start building tst for current fasta
            build_tree_cpp(
                reader.to_string(),
                pam.pamseq,
                pam.guide_size + pam.size,
                pam.size,
                pam.upstream,
                threads,
            )
        except Exception as e:
            exception_handler(
                CrispritzTstError,
                f"Failed building ternary search tree on {fasta}",
                os.EX_DATAERR,
                debug,
                e,
            )
