""" """

from .crispritz_errors import CrispritzPamError
from .exception_handlers import exception_handler

from typing import Tuple

import os


class PAM:

    def __init__(self, fname: str, debug: bool) -> None:
        self._debug = debug  # store debug flag
        # set pam sequence and size
        pam_seq, size = self._read_pam_file(fname)
        self._upstream = size < 0  # pam is upstream
        self._size = abs(size)  # adjust size
        self._guide_size = len(pam_seq) - self._size
        self._pam_seq = self._refine_pam(pam_seq)

    def __len__(self) -> int:
        return self._size

    def _read_pam_file(self, fname: str) -> Tuple[str, int]:
        # expected PAM file format:
        # NNNNNNNNNNNNNNNNNNNNNGG 3 (negative if pam upstream wrt grna)
        try:
            with open(fname, mode="r") as fin:
                pamseq, size = fin.readline().strip().split()
        except (IOError, Exception) as e:
            exception_handler(
                CrispritzPamError,
                f"Failed parsing PAM file {fname}",
                os.EX_IOERR,
                self._debug,
                e,
            )
        if size == 0:  # invalid pam size
            exception_handler(
                CrispritzPamError,
                f"Invalid PAM size: {size}",
                os.EX_DATAERR,
                self._debug,
            )
        return pamseq, int(size)

    def _refine_pam(self, pam_seq: str) -> str:
        # return just the pam sequence
        return pam_seq[: self._size] if self._upstream else pam_seq[-self._size :]

    @property
    def pamseq(self) -> str:
        return self._pam_seq

    @property
    def upstream(self) -> bool:
        return self._upstream

    @property
    def size(self) -> int:
        return self._size

    @property
    def guide_size(self) -> int:
        return self._guide_size
