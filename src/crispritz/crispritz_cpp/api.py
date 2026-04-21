from crispritz import _ternary_search_tree as tst  # type: ignore


def build_tree_cpp(
    sequence: str,
    contig: str,
    pam: str,
    pam_length: int,
    pam_size: int,
    upstream: bool,
    max_bulges: int = 0,
    threads: int = 1,
) -> None:
    """Call the C++ TST builder for a single chromosome sequence.

    Parameters
    ----------
    sequence:
        Full genomic sequence (single chromosome, uppercase IUPAC).
    contig:
        Chromosome / contig identifier used in the output filename(s).
    pam:
        PAM-only string (e.g. ``"NGG"``), without guide placeholder Ns.
    pam_length:
        Total length of the PAM + guide pattern.
    pam_size:
        Length of the PAM portion only.
    upstream:
        True when the PAM precedes the guide (e.g. Cas12a).
    max_bulges:
        Maximum number of bulges allowed during index construction.
    threads:
        Number of OpenMP threads for the PAM search phase.
    """
    tst.build_tree(
        sequence, contig, pam, pam_length, pam_size, upstream, max_bulges, threads
    )
