""" """

from ..crispritz_errors import CrispritzTstError
from ..exception_handlers import exception_handler
from ..crispritz_cpp import build_tree_cpp
from ..verbosity import VERBOSITY_LVL, print_verbosity
from ..genome_io import GenomeReader
from ..pam import PAM

from typing import List

import os


def build_ternary_search_tree(
    fastas: List[str],
    pam_file: str,
    bmax: int,
    outdir: str,
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
    outdir:
        Path to the directory where the genome index will be stored.
    threads:
        Number of OpenMP threads for the PAM search phase.
    verbosity:
        Verbosity level (see ``VERBOSITY_LVL``).
    debug:
        When *True*, exceptions propagate with full stack traces.
    """
    pam = PAM(pam_file, debug)
    for fasta in fastas:
        reader = GenomeReader(fasta, debug)
        reader.read()
        # the contig name is used in the output .bin filename(s).
        # Add leading 'chr' prefix to improve the legacy naming used by
        # (even though the search binary expects e.g. "1" not "chr1").
        contig = (
            reader.header
            if reader.header.startswith("chr")
            else f"chr{reader.header}"
        )
        print_verbosity(
            f"Building TST index for {contig} ({fasta})", verbosity, VERBOSITY_LVL[2]
        )
        try:
            build_tree_cpp(
                reader.to_string(),
                contig,
                pam.pamseq,
                pam.guide_size + pam.size,
                pam.size,
                pam.upstream,
                outdir,
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
