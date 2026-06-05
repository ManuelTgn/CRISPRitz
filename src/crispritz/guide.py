""" """

from .crispritz_errors import CrispritzGuideError
from .dna_alphabet import DNA
from .exception_handlers import exception_handler
from .pam import PAM
from .dna_alphabet import reverse_complement

from typing import List

import os


class Guide:

    def __init__(self, sequence: str, pam: PAM, debug: bool) -> None:
        self._debug = debug  # store debug flag
        # set guide sequence and its reverse complement
        self._sequence = _validate_guide_sequence(
            sequence, len(pam), pam.upstream, self._debug
        )
        self._sequence_rc = reverse_complement(self._sequence)

    def __len__(self) -> int:
        return len(self._sequence)

    @property
    def sequence(self) -> str:
        return self._sequence
    
    @property
    def reverse(self) -> str:
        return self._sequence_rc


class GuideList:
    def __init__(self, fname: str, pam: PAM, debug: bool) -> None:
        self._debug = debug  # store debug flag
        # read guides in input file
        self._guides = self._read_guides_file(pam, fname)

    def _read_guides_file(self, pam: PAM, fname: str) -> List[Guide]:
        # expected guides file format (2 guides):
        # CTAACAGTTGCTTTTATCACNNN
        # CTAACAGTTGCTTTTATCACNNN
        try:
            with open(fname, mode="r") as fin:
                guides = [Guide(line.strip(), pam, self._debug) for line in fin]
        except (IOError, Exception) as e:
            exception_handler(
                CrispritzGuideError,
                f"Failed parsing guide file {fname}",
                os.EX_IOERR,
                self._debug,
                e,
            )
        if not guides:  # no guide found in input file?
            exception_handler(
                CrispritzGuideError,
                f"No guide found in guides file {fname}",
                os.EX_IOERR,
                self._debug,
            )
        # guides must have same length
        if any(len(g) != guides[0] for g in guides[1:]):
            exception_handler(
                CrispritzGuideError,
                "Mismarching guide lengths",
                os.EX_DATAERR,
                self._debug,
            )
        return guides

    @property
    def guides(self) -> List[Guide]:
        return self._guides


def _validate_guide_sequence(
    guide: str, pam_size: int, pam_upstream: bool, debug: bool
) -> str:
    # extract pam slice from guide sequence: should all be 'N'
    expected_pam_slice = guide[:pam_size] if pam_upstream else guide[-pam_size:]
    if any(nt != DNA[4] for nt in expected_pam_slice):
        exception_handler(
            CrispritzGuideError,
            f"Bad guide formatting, guide: {guide}",
            os.EX_DATAERR,
            debug,
        )
    return guide
