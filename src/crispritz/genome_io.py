""" """

from .crispritz_errors import GenomeReaderError, GenomeWriterError
from .exception_handlers import exception_handler
from .enrichment.variants import IndelPair

from typing import List

import os


# define indel upstream and downstream offset
INDELOFFSET = 50


class GenomeReader:

    def __init__(self, fasta_path: str, debug: bool):
        self._debug = debug  # store debug flag
        self._fasta_path = fasta_path
        self._sequence: List[str] = []
        self._sequence_enr: List[str] = []
        self._header: str = ""

    def __repr__(self) -> str:
        seq_preview = f"[{len(self._sequence)} bases]" if self._sequence else "None"
        return (
            f"GenomeReader(fasta_path={self._fasta_path!r}, "
            f"header={self._header!r}, "
            f"sequence={seq_preview}, "
            f"debug={self._debug!r})"
        )

    def __str__(self) -> str:
        if not self._header or not self._sequence:
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
        self._sequence_enr = list(self._sequence)

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

    def insert_snp(self, iupac_nt: str, pos: int) -> None:
        if not self._sequence_enr:
            exception_handler(
                GenomeReaderError,
                "Sequence not initialized, impossible to insert SNPs",
                os.EX_DATAERR,
                self._debug,
            )
        self._sequence_enr[pos] = iupac_nt  # add snp as iupac to sequence

    def insert_indel(self, indel: str, pos: int, offset: int) -> IndelPair:
        refseq = self._sequence[pos - INDELOFFSET : pos + INDELOFFSET + offset]
        indelseq = refseq[:INDELOFFSET] + list(indel) + refseq[INDELOFFSET + offset :]
        return IndelPair(refseq=refseq, indelseq=indelseq)

    def to_string(self) -> str:
        if not self._sequence_enr:
            exception_handler(
                GenomeReaderError,
                f"Sequence not initialized, impossible to retrieve sequence as {str.__name__}",
                os.EX_DATAERR,
                self._debug,
            )
        return "".join(self._sequence_enr)

    @property
    def header(self) -> str:
        assert self._header  # if used, always initialized
        return self._header

    @property
    def sequence(self) -> List[str]:
        assert self._sequence  # if used, always initialized
        return self._sequence

    @property
    def sequence_enr(self) -> List[str]:
        assert self._sequence_enr  # if used, always initialized
        return self._sequence_enr

    @property
    def fname(self) -> str:
        return self._fasta_path


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
