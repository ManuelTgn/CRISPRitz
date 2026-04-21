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
    """Build a Ternary Search Tree index for every input FASTA file.

    Reads each FASTA, extracts the chromosome name from the header, and
    calls the C++ TST builder.  One or more ``.bin`` partition files are
    written per chromosome into the current working directory.

    Parameters
    ----------
    fastas:
        Paths to the per-chromosome FASTA files to index.
    pam_file:
        Path to the PAM specification file.
    bmax:
        Maximum number of bulges; passed to the C++ builder so it extracts
        enough extra bases per site to support bulge-aware search.
    threads:
        Number of OpenMP threads for the PAM search phase.
    verbosity:
        Verbosity level (see ``VERBOSITYLVL``).
    debug:
        When *True*, exceptions propagate with full stack traces.
    """
    pam = PAM(pam_file, debug)
    for fasta in fastas:
        reader = GenomeReader(fasta, debug)
        reader.read()
        # the contig name is used in the output .bin filename(s).
        # Strip any leading 'chr' prefix to match the legacy naming used by
        # mainParallel.cpp (the search binary expects e.g. "1" not "chr1").
        contig = (
            reader.header.replace("chr", "")
            if reader.header.startswith("chr")
            else reader.header
        )
        print_verbosity(
            f"Building TST index for {contig} ({fasta})", verbosity, VERBOSITYLVL[2]
        )
        try:
            build_tree_cpp(
                reader.to_string(),
                contig,
                pam.pamseq,
                pam.guide_size + pam.size,
                pam.size,
                pam.upstream,
                bmax,
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
