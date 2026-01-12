""" """

from .crispritz_error import GenomeReaderError, GenomeWriterError
from .exception_handlers import exception_handler
from .utils import find_fasta_index

from itertools import groupby
from pysam import faidx, FastaFile
from typing import Tuple, List, Optional
from pathlib import Path

import os


class GenomeReader:

    def __init__(self, fasta_path: str, debug: bool):
        self._debug = debug  # store debug flag
        self._fasta_path = fasta_path
        self._sequence: Optional[List[str]] = None
        self._header: Optional[str] = None

    def __repr__(self) -> str:
        seq_preview = (
            f"[{len(self._sequence)} bases]" if self._sequence is not None else "None"
        )
        return (
            f"GenomeReader(fasta_path={self._fasta_path!r}, "
            f"header={self._header!r}, "
            f"sequence={seq_preview}, "
            f"debug={self._debug!r})"
        )

    def __str__(self) -> str:
        if self._header is None or self._sequence is None:
            status = "not read"
            details = ""
        else:
            seq_len = len(self._sequence)
            seq_preview = "".join(self._sequence[:50])
            if seq_len > 50:
                seq_preview += f"... ({seq_len - 50} more bases)"
            status = "read"
            details = f"\n  Header: {self._header}\n  Length: {seq_len} bases\n  Preview: {seq_preview}"
        return f"GenomeReader: {self._fasta_path} ({status}){details}"

    def _extract_header(self, header: str) -> None:
        self._header = header
        if not self._header:
            exception_handler(
                GenomeReaderError,
                f"FASTA header is empty: {self._fasta_path}",
                os.EX_IOERR,
                self._debug,
            )

    def _extract_sequence(self, lines: List[str]) -> None:
        sequence_lines = [line for line in lines if not line.startswith(">")]
        if not sequence_lines:
            exception_handler(
                GenomeReaderError,
                f"FASTA file contains no sequence data: {self._fasta_path}",
                os.EX_IOERR,
                self._debug,
            )
        # join and convert to list of characters
        self._sequence = list("".join(sequence_lines).upper())

    def read(self) -> None:
        try:
            with open(self._fasta_path, mode="r") as fin:
                lines = [line.rstrip("\n\r") for line in fin if line.strip()]
                if not lines:
                    exception_handler(
                        GenomeReaderError,
                        f"FASTA file is empty: {self._fasta_path}",
                        os.EX_IOERR,
                        self._debug,
                    )
                if not lines[0].startswith(">"):
                    exception_handler(
                        GenomeReaderError,
                        f"FASTA file must start with a header line: {self._fasta_path}",
                        os.EX_IOERR,
                        self._debug,
                    )
                # extract header (remove '>' and strip whitespace)
                self._extract_header(lines[0].lstrip(">").strip())
                # extract sequence (all non-header lines)
                self._extract_sequence(lines[1:])
        except (IOError, Exception) as e:
            exception_handler(
                GenomeReaderError,
                f"Failed reading FASTA: {self._fasta_path}",
                os.EX_IOERR,
                self._debug,
                e,
            )

    @property
    def header(self) -> Optional[str]:
        return self._header

    @property
    def sequence(self) -> Optional[List[str]]:
        return self._sequence


class GenomeWriter:

    def __init__(self, outfile: str, debug: bool) -> None:
        self._debug = debug  # store debug flag
        self._outfile = outfile

    def __repr__(self) -> str:
        return f"GenomeWriter(outfile={self._outfile!r}, debug={self._debug!r})"

    def __str__(self) -> str:
        return f"GenomeWriter: {self._outfile}"

    def write(self, header: str, sequence_list: List[str]) -> None:
        sequence = "".join(sequence_list)
        try:
            with open(self._outfile, mode="w") as fout:
                fout.write(f">{header}\n")  # write header
                fout.write(f"{sequence}\n")  # write sequence with newline
        except (IOError, Exception) as e:
            exception_handler(
                GenomeWriterError,
                f"Failed writing FASTA: {self._outfile}",
                os.EX_IOERR,
                self._debug,
                e,
            )
