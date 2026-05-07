from crispritz import _ternary_search_tree as tst  # type: ignore


def build_tree_cpp(
    sequence: str,
    contig: str,
    pam: str,
    pam_length: int,
    pam_size: int,
    upstream: bool,
    outdir: str,
    max_bulges: int = 0,
    threads: int = 1,
) -> None:
    """Call the C++ TST builder for a single chromosome sequence.

    Parameters
    ----------
    sequence: str
        Full genomic sequence (single chromosome, uppercase IUPAC).
    contig: str
        Chromosome / contig identifier used in the output filename(s).
    pam: str
        PAM-only string (e.g. ``"NGG"``), without guide placeholder Ns.
    pam_length: int
        Total length of the PAM + guide pattern.
    pam_size: int
        Length of the PAM portion only.
    upstream: bool
        True when the PAM precedes the guide (e.g. Cas12a).
    outdir: str
        Path to the directory where the genome index will be stored.
    max_bulges: int
        Maximum number of bulges allowed during index construction.
    threads: int
        Number of OpenMP threads for the PAM search phase.
    """
    tst.build_tree(
        sequence, contig, pam, pam_length, pam_size, upstream, outdir, max_bulges, threads
    )
