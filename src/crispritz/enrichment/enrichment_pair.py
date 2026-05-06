""" """

from .crispritz_enrichment_error import EnrichmentPairError

from ..exception_handlers import exception_handler

import os


class EnrichPair:

    def __init__(self, debug: bool) -> None:
        self._debug = debug  # store debug flag
        self._fasta: str = ""
        self._vcf: str = ""

    def __repr__(self) -> str:
        """Return a detailed string representation for debugging."""
        return f"EnrichPair(fasta={self._fasta!r}, vcf={self._vcf!r})"

    def __str__(self) -> str:
        """Return a human-readable string representation."""
        fasta_str = self._fasta if self._fasta is not None else "not set"
        vcf_str = self._vcf if self._vcf is not None else "not set"
        return f"EnrichPair: FASTA={fasta_str}, VCF={vcf_str}"

    @property
    def fasta(self) -> str:
        return self._fasta

    @fasta.setter
    def fasta(self, value: str) -> None:
        if not isinstance(value, str) or not value:
            exception_handler(
                EnrichmentPairError,
                "FASTA file must be a non-empty str, got "
                f"{type(value).__name__} instead",
                os.EX_DATAERR,
                self._debug,
            )
        self._fasta = value

    @property
    def vcf(self) -> str:
        return self._vcf

    @vcf.setter
    def vcf(self, value: str) -> None:
        if not isinstance(value, str):
            exception_handler(
                EnrichmentPairError,
                "VCF file must be a non-empty str, got "
                f"{type(value).__name__} instead",
                os.EX_DATAERR,
                self._debug,
            )
        self._vcf = value
