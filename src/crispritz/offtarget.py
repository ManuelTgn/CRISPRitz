"""Representation of a single CRISPR off-target site.

This module defines :class:`OffTarget`, the canonical in-memory representation
of one off-target candidate produced by the TST search.  Every instance maps
directly to one row of the final TSV report and carries all fields needed for
downstream scoring (CFD, CRISTA, ...).
"""

from __future__ import annotations

from typing import List, Optional, Tuple


# ==============================================================================
# Constants
# ==============================================================================

#: Sentinel used for scores that have not yet been computed.
_SCORE_NA: str = "NA"

#: TSV column order – used by :py:meth:`OffTarget.to_tsv_row`.
TSV_HEADER: List[str] = [
    "chrom",
    "pos",
    "strand",
    "grna",
    "spacer",
    "mismatches",
    "bulge_dna",
    "bulge_rna",
    "bulge_type",
    "cfd_score",
]


# ==============================================================================
# OffTarget class
# ==============================================================================


class OffTarget:
    """Represent a single CRISPR off-target site and its associated scores.

    Each instance corresponds to one candidate off-target locus produced by the
    TST near-search.  The object is intentionally lightweight: it stores raw
    alignment data and lazily-populated scores, and knows how to serialise
    itself to a TSV row for the final report.

    Attributes
    ----------
    _chrom : str
        Chromosome (or contig) name of the off-target site.
    _pos : int
        1-based genomic start position.
    _strand : str
        Strand on which the site was found: ``'+'`` or ``'-'``.
    _grna : str
        Aligned guide RNA sequence including PAM placeholder ``N``.
        Mismatch positions are represented in upper-case; gap characters
        (``'-'``) mark bulge positions.
    _spacer : str
        Aligned genomic target (spacer) sequence including PAM bases.
        Mismatch bases are lower-case (convention from the C++ search).
    _mm : int
        Number of substitution mismatches between guide and target.
    _bulge_type : str
        Bulge classification: ``'X'`` (none), ``'DNA'``, ``'RNA'``, or
        ``'DNA,RNA'``.
    _bdna : int
        Number of DNA-bulge bases (gap in guide / extra base in target).
    _brna : int
        Number of RNA-bulge bases (gap in target / extra base in guide).
    _cfd_score : str
        Cutting Frequency Determination score, or ``'NA'`` when not yet
        computed.  Stored as a string to allow ``'NA'`` as a sentinel
        without requiring ``Optional[float]`` throughout the pipeline.
    _debug : bool
        When ``True``, exceptions propagate with full tracebacks instead of
        being swallowed into user-friendly error messages.

    Examples
    --------
    >>> ot = OffTarget(
    ...     chrom="chr1",
    ...     pos=123456,
    ...     strand="+",
    ...     grna="ACGTACGTACGTACGTACGTNGG",
    ...     spacer="ACGTACGTACGTACGTACGTaGG",
    ...     mm=1,
    ...     bulge_type="X",
    ...     bdna=0,
    ...     brna=0,
    ... )
    >>> ot.cfd_score = "0.85"
    >>> print(ot.to_tsv_row())
    chr1\t123456\t+\tACGTACGTACGTACGTACGTNGG\tACGTACGTACGTACGTACGTaGG\t1\tX\t0\t0\t0.85
    """

    # ==========================================================================
    # Construction
    # ==========================================================================

    def __init__(
        self,
        chrom: str,
        pos: int,
        strand: str,
        grna: str,
        spacer: str,
        mm: int,
        bulge_type: str,
        bdna: int,
        brna: int,
        cfd_score: str = _SCORE_NA,
        debug: bool = False,
    ) -> None:
        """Initialise an :class:`OffTarget` instance.

        Parameters
        ----------
        chrom:
            Chromosome name (e.g. ``'chr1'``).
        pos:
            1-based genomic start position (must be > 0).
        strand:
            Strand indicator; must be ``'+'`` or ``'-'``.
        grna:
            Aligned guide RNA sequence with PAM (``N``-placeholder) appended
            or prepended depending on the PAM orientation.
        spacer:
            Aligned genomic target sequence with PAM bases.  Lower-case
            characters indicate mismatch positions; ``'-'`` indicates bulge
            gaps.
        mm:
            Number of substitution mismatches (>= 0).
        bulge_type:
            One of ``'X'``, ``'DNA'``, ``'RNA'``, ``'DNA,RNA'``.
        bdna:
            Number of DNA-bulge bases (>= 0).
        brna:
            Number of RNA-bulge bases (>= 0).
        cfd_score:
            Pre-computed CFD score string, or ``'NA'`` (default).
        debug:
            Enable debug mode for verbose error propagation.

        Raises
        ------
        ValueError
            If any argument fails basic sanity checks (invalid strand, negative
            counts, unknown bulge type, …).
        """
        self._debug = debug
        self._chrom = chrom
        self._pos = pos
        self._strand = strand
        self._grna = grna
        self._spacer = spacer
        self._mm = mm
        self._bulge_type = bulge_type
        self._bdna = bdna
        self._brna = brna
        self._cfd_score: str = cfd_score

    # ------------------------------------------------------------------
    # Representation
    # ------------------------------------------------------------------

    def __repr__(self) -> str:  # noqa: D105
        return (
            f"OffTarget("
            f"chrom={self._chrom!r}, "
            f"pos={self._pos}, "
            f"strand={self._strand!r}, "
            f"mm={self._mm}, "
            f"bulge_type={self._bulge_type!r}, "
            f"bdna={self._bdna}, "
            f"brna={self._brna}, "
            f"cfd={self._cfd_score!r}"
            f")"
        )

    def __str__(self) -> str:  # noqa: D105
        return (
            f"{self._chrom}:{self._pos}{self._strand}  "
            f"gRNA={self._grna}  spacer={self._spacer}  "
            f"mm={self._mm}  bulge={self._bulge_type}({self._bdna}/{self._brna})  "
            f"CFD={self._cfd_score}"
        )

    def __eq__(self, other: object) -> bool:
        """Two off-targets are equal when they map to the same genomic locus.

        Identity is defined by chromosome, position, strand, and the aligned
        spacer sequence (which fully determines the target window).  Score
        fields are intentionally excluded so that the same hit scored by
        different methods compares equal.

        Parameters
        ----------
        other:
            Object to compare against.

        Returns
        -------
        bool
            ``True`` when *other* is an :class:`OffTarget` with the same locus.
        """
        if not isinstance(other, OffTarget):
            return NotImplemented
        return (
            self._chrom == other._chrom
            and self._pos == other._pos
            and self._strand == other._strand
            and self._spacer == other._spacer
        )

    def __hash__(self) -> int:
        """Hash based on the locus identity (mirrors :py:meth:`__eq__`).

        Returns
        -------
        int
            Hash value suitable for use in sets and dict keys.
        """
        return hash((self._chrom, self._pos, self._strand, self._spacer))

    # ==========================================================================
    # Properties – getters
    # ==========================================================================

    @property
    def chrom(self) -> str:
        """Chromosome name of the off-target site."""
        return self._chrom

    @property
    def pos(self) -> int:
        """1-based genomic start position."""
        return self._pos

    @property
    def strand(self) -> str:
        """Strand of the off-target site (``'+'`` or ``'-'``)."""
        return self._strand

    @property
    def grna(self) -> str:
        """Aligned guide RNA sequence (with PAM placeholder Ns)."""
        return self._grna

    @property
    def spacer(self) -> str:
        """Aligned genomic target (spacer) sequence including PAM bases."""
        return self._spacer

    @property
    def mm(self) -> int:
        """Number of substitution mismatches."""
        return self._mm

    @property
    def bulge_type(self) -> str:
        """Bulge classification (``'X'``, ``'DNA'``, ``'RNA'``, or ``'DNA/RNA'``)."""
        return self._bulge_type

    @property
    def bdna(self) -> int:
        """Number of DNA-bulge bases (gap in guide)."""
        return self._bdna

    @property
    def brna(self) -> int:
        """Number of RNA-bulge bases (gap in target)."""
        return self._brna

    @property
    def total_score(self) -> int:
        """Total edit distance: mismatches + DNA bulges + RNA bulges.

        Returns
        -------
        int
            ``mm + bdna + brna``.
        """
        return self._mm + self._bdna + self._brna

    @property
    def has_bulge(self) -> bool:
        """``True`` when the off-target involves at least one bulge position.

        Returns
        -------
        bool
            ``True`` if *bulge_type* is not ``'X'``.
        """
        return self._bulge_type != "X"

    @property
    def cfd_score(self) -> str:
        """CFD score string, or ``'NA'`` when not yet computed."""
        return self._cfd_score

    @cfd_score.setter
    def cfd_score(self, value: float) -> None:
        """Set the CFD score.

        Parameters
        ----------
        value:
            A numeric string (e.g. ``'0.85'``) or ``'NA'``.

        Raises
        ------
        TypeError
            If *value* is not a :class:`str`.
        """
        if not isinstance(value, float):
            raise TypeError(
                f"cfd_score must be a str, got {type(value).__name__!r}"
            )
        self._cfd_score = f"{value:.2f}"

    @property
    def debug(self) -> bool:
        """Debug mode flag."""
        return self._debug

    # ==========================================================================
    # Serialisation
    # ==========================================================================

    def to_tsv_row(self, sep: str = "\t") -> str:
        """Serialise the off-target to a single delimited row.

        The column order matches :data:`TSV_HEADER` and is stable across all
        :class:`OffTarget` instances, making it safe to concatenate rows
        directly without a per-row header.

        Parameters
        ----------
        sep:
            Field separator (default: tab ``'\\t'``).

        Returns
        -------
        str
            A single line **without** a trailing newline.

        Examples
        --------
        >>> row = ot.to_tsv_row()
        >>> fields = row.split("\\t")
        >>> assert len(fields) == len(TSV_HEADER)
        """
        fields = [
            self._chrom,
            str(self._pos),
            self._strand,
            self._grna,
            self._spacer,
            str(self._mm),
            self._bulge_type,
            str(self._bdna),
            str(self._brna),
            self._cfd_score,
        ]
        return sep.join(fields)

    # ==========================================================================
    # Convenience helpers
    # ==========================================================================

    def cfd_score_as_float(self) -> Optional[float]:
        """Return the CFD score as a :class:`float`, or ``None`` if not set.

        Returns
        -------
        float or None
            Parsed CFD value, or ``None`` when :py:attr:`cfd_score` is
            ``'NA'`` or cannot be parsed as a float.
        """
        if self._cfd_score == _SCORE_NA:
            return None
        try:
            return float(self._cfd_score)
        except ValueError:
            return None

    def locus(self) -> str:
        """Return a compact human-readable locus string.

        Returns
        -------
        str
            E.g. ``'chr1:123456(+)'``.
        """
        return f"{self._chrom}:{self._pos}({self._strand})"
    